#! /usr/bin/env python

import os
import numpy as np
import time
from datetime import date
import logging
import csv
import networkx as nx
from node2vec import Node2Vec
from cellmaps_utils import constants
from cellmaps_utils import logutils
from cellmaps_utils.provenance import ProvenanceUtil
import warnings

import cellmaps_ppi_embedding
from cellmaps_ppi_embedding.exceptions import CellMapsPPIEmbeddingError

logger = logging.getLogger(__name__)


class EmbeddingGenerator(object):
    """
    Base class for implementations that generate
    network embeddings
    """
    DIMENSIONS = 1024

    def __init__(self, dimensions=DIMENSIONS):
        """
        Constructor
        """
        self._dimensions = dimensions

    def get_dimensions(self):
        """
        Gets number of dimensions this embedding will generate

        :return: number of dimensions aka vector length
        :rtype: int
        """
        return self._dimensions

    def get_next_embedding(self):
        """
        Generator method for getting next embedding.
        Caller should implement with ``yield`` operator

        :raises: NotImplementedError: Subclasses should implement this
        :return: Embedding
        :rtype: list
        """
        raise NotImplementedError('Subclasses should implement')


class Node2VecEmbeddingGenerator(EmbeddingGenerator):
    """

    """
    P_DEFAULT = 2
    Q_DEFAULT = 1
    WALK_LENGTH = 80
    NUM_WALKS = 10
    WORKERS = 8

    def __init__(self, nx_network, p=P_DEFAULT, q=Q_DEFAULT, dimensions=EmbeddingGenerator.DIMENSIONS,
                 walk_length=WALK_LENGTH, num_walks=NUM_WALKS, workers=WORKERS):
        """
        Constructor
        """
        super().__init__(dimensions=dimensions)
        self._nx_network = nx_network
        self._p = p
        self._q = q
        self._walk_length = walk_length
        self._num_walks = num_walks
        self._workers = workers

    def _remove_header_edge_from_network(self):
        """
        Removes edge named geneA geneB which is actually the
        header of the edge list file

        """
        # remove geneA geneB edge cause it is the header of file
        try:
            self._nx_network.remove_edge('geneA', 'geneB')
        except nx.NetworkXError as ne:
            logger.debug('No edge named geneA -> geneB to remove' + str(ne))

        self._nx_network.remove_nodes_from(['geneA', 'geneB'])

    def get_next_embedding(self):
        """

        :return:
        """
        if self._nx_network is None:
            raise CellMapsPPIEmbeddingError('network is None')

        self._remove_header_edge_from_network()

        n2v_obj = Node2Vec(self._nx_network, dimensions=self._dimensions,
                           walk_length=self._walk_length,
                           num_walks=self._num_walks,
                           workers=self._workers, q=self._q, p=self._p)

        # Embed nodes
        model = n2v_obj.fit(window=10, min_count=0, sg=1, epochs=1)
        for key in model.wv.index_to_key:
            row = [key]
            row.extend(model.wv[key].tolist())
            yield row


class FakeEmbeddingGenerator(EmbeddingGenerator):
    """
    Fakes PPI embedding
    """

    def __init__(self, ppi_downloaddir, dimensions=1024):
        """
        Constructor

        :param dimensions: Desired size of output embedding
        :type dimensions: int
        """
        super().__init__(dimensions=dimensions)

        self._ppi_downloaddir = ppi_downloaddir
        self._gene_list = self._get_gene_list()

        warnings.warn(constants.PPI_EMBEDDING_FILE +
                      ' contains FAKE DATA!!!!\n'
                      'You have been warned\nHave a nice day\n')
        logger.error(constants.PPI_EMBEDDING_FILE +
                     ' contains FAKE DATA!!!! '
                     'You have been warned. Have a nice day')

    def _get_gene_list(self):

        ppi_gene_node_attrs_file = os.path.join(self._ppi_downloaddir, constants.PPI_GENE_NODE_ATTR_FILE)
        gene_list = []

        with open(ppi_gene_node_attrs_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_list.append(row['name'])
        return gene_list

    def get_next_embedding(self):
        """
        Generator method for getting next embedding.
        Caller should implement with ``yield`` operator

        :raises: NotImplementedError: Subclasses should implement this
        :return: Embedding
        :rtype: list
        """
        for g in self._gene_list:
            row = [g]
            row.extend(np.random.normal(size=self.get_dimensions()))  # sample normal distribution
            yield row


class CellMapsPPIEmbedder(object):
    """
    Executes the generation of protein-protein interaction (PPI) embeddings.

    This class coordinates reading PPI data, generating node embeddings using a
    specified embedding generator (e.g., Node2Vec), writing results to disk, and
    registering datasets and software in FAIRSCAPE.
    """
    PPI_EDGELIST_FILEKEY = 'edgelist'

    def __init__(self, outdir=None,
                 embedding_generator=None,
                 inputdir=None,
                 skip_logging=True,
                 name=None,
                 organization_name=None,
                 project_name=None,
                 provenance_utils=ProvenanceUtil(),
                 input_data_dict=None,
                 provenance=None):
        """
        Constructor

        :param outdir: directory where ppi embeddings will be saved
        :type outdir: str
        :param embedding_generator: Object responsible for generating the embeddings.
                                    Must implement `get_next_embedding()`, typically an
                                    instance of a subclass of :py:class:`~EmbeddingGenerator`, such as
                                    :py:class:`~Node2VecEmbeddingGenerator` or :py:class:`~FakeEmbeddingGenerator`.
        :type embedding_generator: :py:class:`~EmbeddingGenerator`
        :param inputdir: Input directory that contains ppi edgelist file and its RO-Crate metadata file.
        :type inputdir: str or None
        :param skip_logging: If ``True`` skip logging, if ``None`` or ``False`` do NOT skip logging
        :type skip_logging: bool
        :param name: Optional display name for the dataset. If not provided, the name will
                     be inferred from the RO-Crate metadata or provenance dictionary.
        :type name: str or None
        :param organization_name: Optional name of the organization generating the dataset.
                          Used in provenance tracking. Falls back to RO-Crate or provenance input if missing.
        :type organization_name: str or None
        :param project_name: Optional name of the project associated with the dataset.
                     Used in provenance tracking. Falls back to RO-Crate or provenance input if missing.
        :type project_name: str or None
        :param provenance_utils: Utility class used for RO-Crate generation and FAIRSCAPE
                                 dataset, software, and computation registration. Defaults to a new
                                 :py:class:`~ProvenanceUtil`.
        :type provenance_utils: :py:class:`~ProvenanceUtil`
        :param input_data_dict: Dictionary of parameters and their values that capture the
                                configuration used to generate the embeddings. This is serialized
                                in the task metadata for reproducibility. If not provided, one is
                                auto-generated from available parameters.

                                Example:

                                .. code-block:: python

                                    {'outdir': '/output/path', 'inputdir': '/input/path'}
        :type input_data_dict: dict or None
        :param provenance: Optional dictionary specifying provenance metadata. Required if
                           `inputdir` does not contain an RO-Crate. Used to describe the
                           input edgelist, dataset authorship, and context.

                           Example:

                           .. code-block:: python

                               {
                                   'name': 'Example PPI Dataset',
                                   'organization-name': 'CM4AI',
                                   'project-name': 'Network Embedding',
                                   'description': 'Node2Vec embeddings of protein-protein interactions',
                                   'keywords': ['PPI', 'embedding', 'node2vec'],
                                   'edgelist': {
                                       'name': 'PPI Edgelist File',
                                       'author': 'Krogan Lab',
                                       'version': '1.0',
                                       'data-format': 'tsv'
                                   }
                               }
        :type provenance: dict or None
        """
        if outdir is None:
            raise CellMapsPPIEmbeddingError('outdir is None')

        self._outdir = os.path.abspath(outdir)
        self._inputdir = os.path.abspath(inputdir) if inputdir is not None else inputdir
        self._start_time = int(time.time())
        self._end_time = -1
        self._embedding_generator = embedding_generator
        self._name = name
        self._project_name = project_name
        self._organization_name = organization_name
        self._keywords = None
        self._description = None
        self._input_data_dict = input_data_dict
        self._provenance_utils = provenance_utils
        self._provenance = provenance
        self._inputdataset_ids = []
        if skip_logging is None:
            self._skip_logging = False
        else:
            self._skip_logging = skip_logging

        if self._input_data_dict is None or not self._input_data_dict:
            self._input_data_dict = {'outdir': self._outdir,
                                     'inputdir': self._inputdir,
                                     'embedding_generator': str(self._embedding_generator),
                                     'name': self._name,
                                     'project_name': self._project_name,
                                     'organization_name': self._organization_name,
                                     'skip_logging': self._skip_logging,
                                     'provenance': str(self._provenance)
                                     }

        logger.debug('In constructor')

    @staticmethod
    def get_apms_edgelist_file(input_dir=None,
                               edgelist_filename=constants.PPI_EDGELIST_FILE):
        """
        :param edgelist_filename:
        :return:
        :param input_dir:
        :return:
        """
        return os.path.join(input_dir, edgelist_filename)

    def _write_task_start_json(self):
        """
        Writes task_start.json file with information about
        what is to be run

        """
        logutils.write_task_start_json(outdir=self._outdir,
                                       start_time=self._start_time,
                                       version=cellmaps_ppi_embedding.__version__,
                                       data={'commandlineargs': self._input_data_dict})

    def _update_provenance_fields(self):
        """

        :return:
        """
        if os.path.exists(os.path.join(self._inputdir, constants.RO_CRATE_METADATA_FILE)):
            prov_attrs = self._provenance_utils.get_merged_rocrate_provenance_attrs(self._inputdir,
                                                                                    override_name=self._name,
                                                                                    override_project_name=self._project_name,
                                                                                    override_organization_name=self._organization_name,
                                                                                    extra_keywords=['AP-MS Embedding',
                                                                                                    'AP-MS',
                                                                                                    'embedding'])

            self._name = prov_attrs.get_name()
            self._organization_name = prov_attrs.get_organization_name()
            self._project_name = prov_attrs.get_project_name()
            self._keywords = prov_attrs.get_keywords()
            self._description = prov_attrs.get_description()
        elif self._provenance is not None:
            self._name = self._provenance['name'] if 'name' in self._provenance else 'PPI Embedding'
            self._organization_name = self._provenance['organization-name'] \
                if 'organization-name' in self._provenance else 'NA'
            self._project_name = self._provenance['project-name']\
                if 'project-name' in self._provenance else 'NA'
            self._keywords = self._provenance['keywords'] if 'keywords' in self._provenance else ['ppi']
            self._description = self._provenance['description'] if 'description' in self._provenance else \
                'Embedding of PPIs'
        else:
            raise CellMapsPPIEmbeddingError('Input directory should be an RO-Crate or provenance should be specified.')

    def _create_run_crate(self):
        """
        Creates rocrate for output directory

        :raises CellMapsProvenanceError: If there is an error
        """
        logger.debug('Registering rocrate with FAIRSCAPE')

        try:
            self._provenance_utils.register_rocrate(self._outdir,
                                                    name=self._name,
                                                    organization_name=self._organization_name,
                                                    project_name=self._project_name,
                                                    description=self._description,
                                                    keywords=self._keywords)
        except TypeError as te:
            raise CellMapsPPIEmbeddingError('Invalid provenance: ' + str(te))
        except KeyError as ke:
            raise CellMapsPPIEmbeddingError('Key missing in provenance: ' + str(ke))

    def _register_software(self):
        """
        Registers this tool

        :raises CellMapsImageEmbeddingError: If fairscape call fails
        """
        logger.debug('Registering software with FAIRSCAPE')
        software_keywords = self._keywords
        software_keywords.extend(['tools', cellmaps_ppi_embedding.__name__])
        software_description = self._description + ' ' + \
                               cellmaps_ppi_embedding.__description__
        self._softwareid = self._provenance_utils.register_software(self._outdir,
                                                                    name=cellmaps_ppi_embedding.__name__,
                                                                    description=software_description,
                                                                    author=cellmaps_ppi_embedding.__author__,
                                                                    version=cellmaps_ppi_embedding.__version__,
                                                                    file_format='py',
                                                                    keywords=software_keywords,
                                                                    url=cellmaps_ppi_embedding.__repo_url__)

    def _register_computation(self):
        """
        # Todo: added in used dataset, software and what is being generated
        :return:
        """
        logger.debug('Getting id of input rocrate')
        if os.path.exists(os.path.join(self._inputdir, constants.RO_CRATE_METADATA_FILE)):
            self._inputdataset_ids.append(self._provenance_utils.get_id_of_rocrate(self._inputdir))
        else:
            self._register_input_datasets()

        logger.debug('Registering computation with FAIRSCAPE')
        keywords = self._keywords
        keywords.extend(['computation'])
        description = self._description + ' run of ' + cellmaps_ppi_embedding.__name__
        self._provenance_utils.register_computation(self._outdir,
                                                    name=cellmaps_ppi_embedding.__computation_name__,
                                                    run_by=str(self._provenance_utils.get_login()),
                                                    command=str(self._input_data_dict),
                                                    description=description,
                                                    keywords=keywords,
                                                    used_software=[self._softwareid],
                                                    used_dataset=self._inputdataset_ids,
                                                    generated=[self._embedding_file_id])

    def _register_input_datasets(self):
        """
        Registers cm4ai/apms dataset or samples and unique input
        datasets with FAIRSCAPE
        adding values to **self._inputdataset_ids**

        """
        if CellMapsPPIEmbedder.PPI_EDGELIST_FILEKEY not in self._provenance:
            logger.debug('Provenance information for input file (PPI edgelist) was not specified.')
            return
        ppi_edgelist_datasetid = None
        if 'guid' in self._provenance[CellMapsPPIEmbedder.PPI_EDGELIST_FILEKEY]:
            ppi_edgelist_datasetid = self._provenance[CellMapsPPIEmbedder.PPI_EDGELIST_FILEKEY]['guid']

        if ppi_edgelist_datasetid is not None:
            self._inputdataset_ids.append(ppi_edgelist_datasetid)
            logger.debug('PPI edgelist have dataset id.')
        else:
            ppi_edgelist_datasetid = self._provenance_utils.register_dataset(
                self._outdir, source_file=self.get_apms_edgelist_file(self._inputdir),
                data_dict=self._provenance[CellMapsPPIEmbedder.PPI_EDGELIST_FILEKEY], skip_copy=False)
            self._inputdataset_ids.append(ppi_edgelist_datasetid)
            logger.debug('PPI edgelist dataset id: ' + str(ppi_edgelist_datasetid))

    def _register_embedding_file(self):
        """
        Registers embedding file as a dataset

        :return: id of datafile dataset
        :rtype: str
        """
        logger.debug('Registering embedding file with FAIRSCAPE')
        description = self._description
        description += ' file'
        keywords = self._keywords
        keywords.extend(['file'])
        data_dict = {'name': cellmaps_ppi_embedding.__name__ + ' output file',
                     'description': description,
                     'keywords': keywords,
                     'data-format': 'tsv',
                     'author': cellmaps_ppi_embedding.__name__,
                     'version': cellmaps_ppi_embedding.__version__,
                     'schema': 'https://raw.githubusercontent.com/fairscape/cm4ai-schemas/main/v0.1.0/cm4ai_schema_apms_embedding.json',
                     'date-published': date.today().strftime(self._provenance_utils.get_default_date_format_str())}
        self._embedding_file_id = self._provenance_utils.register_dataset(self._outdir,
                                                                          source_file=self.get_ppi_embedding_file(),
                                                                          data_dict=data_dict)

    def get_ppi_embedding_file(self):
        """
        Gets PPI embedding file in output directory

        :return:
        :rtype: str
        """
        return os.path.join(self._outdir, constants.PPI_EMBEDDING_FILE)

    def generate_readme(self):
        description = getattr(cellmaps_ppi_embedding, '__description__', 'No description provided.')
        version = getattr(cellmaps_ppi_embedding, '__version__', '0.0.0')

        with open(os.path.join(os.path.dirname(__file__), 'readme_outputs.txt'), 'r') as f:
            readme_outputs = f.read()

        readme = readme_outputs.format(DESCRIPTION=description, VERSION=version)
        with open(os.path.join(self._outdir, 'README.txt'), 'w') as f:
            f.write(readme)

    def run(self):
        """
        Run node2vec to create embeddings


        :return:
        """
        logger.debug('In run method')
        exitcode = 99
        try:
            if not os.path.isdir(self._outdir):
                os.makedirs(self._outdir, mode=0o755)

            if self._skip_logging is False:
                logutils.setup_filelogger(outdir=self._outdir,
                                          handlerprefix='cellmaps_ppi_embedding')
            self._write_task_start_json()

            self.generate_readme()

            self._update_provenance_fields()

            self._create_run_crate()

            self._register_software()

            with open(self.get_ppi_embedding_file(), 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                header_line = ['id']
                header_line.extend([x for x in range(self._embedding_generator.get_dimensions())])
                writer.writerow(header_line)
                for row in self._embedding_generator.get_next_embedding():
                    writer.writerow(row)

            self._register_embedding_file()
            self._register_computation()

            exitcode = 0

        finally:
            self._end_time = int(time.time())
            # write a task finish file
            logutils.write_task_finish_json(outdir=self._outdir,
                                            start_time=self._start_time,
                                            end_time=self._end_time,
                                            status=exitcode)
        return exitcode
