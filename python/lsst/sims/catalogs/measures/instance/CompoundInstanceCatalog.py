from __future__ import with_statement
import numpy
from collections import OrderedDict
from lsst.sims.catalogs.generation.db import CompoundCatalogDBObject


class CompoundInstanceCatalog(object):
    """
    This is essentially an InstanceCatalog class meant to wrap together
    several disparate InstanceCatalog instantiations that will ultimately
    be written to the same output catalog.

    You pass the constructor a list of InstanceCatalog classes, a list of
    CatalogDBObject classes, and ObservationMetaData, and an optional SQL constraint.

    The write_catalog method then writes all of the InstanceCatalogs to one
    ASCII file using the same API as InstanceCatalog.write_catalog.
    """

    def __init__(self, instanceCatalogClassList, catalogDBObjectClassList,
                 obs_metadata=None, constraint=None, compoundDBclass = None):
        """
        @param [in] instanceCatalogClassList is a list of the InstanceCatalog
        classes to be combined into one output catalog.

        @param [in] catalogDBObjectClassList is a list of the CatalogDBObject
        classes to be associated with the InstanceCatalog classes in
        instanceCatalogClassList.  There should be one CatalogDBObject class
        for each InstanceCatalogClass.

        @param [in] obs_metadata is the ObservationMetaData describing
        the telescope pointing

        @param [in] constraint is an optional SQL constraint to be applied
        to the database query

        @param [in] compoundDBclass is an optional argument specifying what
        CompoundCatalogDBobject class(es) to use to combine InstanceCatalogs
        that query the same table.  This can be either a single
        ComboundCatalogDBObject class, or a list of classes.  The
        CompoundInstanceCatalog will figure out which InstanceCatalog(s) go with
        which CompoundCatalogDBObject class.  If no CompoundCatalogDBObject class
        corresponds to a given group of InstanceCatalogs, then the base
        CompoundCatalogDBObject class will be used.

        Note: compoundDBclass should be a CompoundCatalogDBObject class.
        Not an instantiation of a CompoundCatalogDBObject class.
        """

        self._compoundDBclass = compoundDBclass
        self._obs_metadata = obs_metadata
        self._dbo_list = catalogDBObjectClassList
        self._ic_list = instanceCatalogClassList
        self._constraint = constraint
        self._active_connections = []

        assigned = [False]*len(self._dbo_list)
        self._dbObjectGroupList = []

        for ix in range(len(self._dbo_list)):
            for row in self._dbObjectGroupList:
                if self.areDBObjectsTheSame(self._dbo_list[ix], self._dbo_list[row[0]]):
                    row.append(ix)
                    assigned[ix] = True
                    break

            if not assigned[ix]:
                new_row = [ix]
                for iy in range(ix):
                    if not assigned[iy]:
                        if self.areDBObjectsTheSame(self._dbo_list[ix], self._dbo_list[iy]):
                            new_row.append(iy)

                self._dbObjectGroupList.append(new_row)

    def areDBObjectsTheSame(self, db1, db2):
        """
        @param [in] db1 is a CatalogDBObject instantiation

        @param [in] db2 is a CatalogDBObject instantiation

        @param [out] a boolean stating whether or not db1 and db2
        query the same table of the same database
        """

        if hasattr(db1, 'host'):
            host1 = db1.host
        else:
            host1 = None

        if hasattr(db2, 'host'):
            host2 = db2.host
        else:
            host2 = None

        if hasattr(db1, 'port'):
            port1 = db1.port
        else:
            port1 = None

        if hasattr(db2, 'port'):
            port2 = db2.port
        else:
            port2 = None

        if hasattr(db1, 'driver'):
            driver1 = db1.driver
        else:
            driver1 = None

        if hasattr(db2, 'driver'):
            driver2 = db2.driver
        else:
            driver2 = None

        if db1.tableid != db2.tableid:
            return False
        if host1 != host2:
            return False
        if db1.database != db2.database:
            return False
        if port1 != port2:
            return False
        if driver1 != driver2:
            return False
        return True

    def find_a_connection(self, dboClass):
        """
        Find an active database connection for a DBObject

        @param [in] dbo is a DBObject class that needs to be connected

        @param [out] returns a connection of self._active_connections
        that suits the DBObject.  Returns None otherwise.
        """

        if hasattr(dboClass, 'database'):
            desired_database = dboClass.database
        else:
            desired_database = None

        if hasattr(dboClass, 'driver'):
            desired_driver = dboClass.driver
        else:
            desired_driver = None

        if hasattr(dboClass, 'host'):
            desired_host = dboClass.host
        else:
            desired_host = None

        if hasattr(dboClass, 'port'):
            desired_port = dboClass.port
        else:
            desired_port = None

        if hasattr(dboClass, 'verbose'):
            desired_verbose = dboClass.verbose
        else:
            desired_verbose = False

        best_connection = None

        for cc in self._active_connections:
            if cc.database is desired_database and \
               cc.driver is desired_driver and \
               cc.host is desired_host and \
               cc.port is desired_port and \
               cc.verbose is desired_verbose:

                best_connection = cc
                break

        return best_connection

    def write_catalog(self, filename, chunk_size=None, write_header=True, write_mode='w'):
        """
        Write the stored list of InstanceCatalogs to a single ASCII output catalog.

        @param [in] filename is the name of the file to be written

        @param [in] chunk_size is an optional parameter telling the CompoundInstanceCatalog
        to query the database in manageable chunks (in case returning the whole catalog
        takes too much memory)

        @param [in] write_header a boolean specifying whether or not to add a header
        to the output catalog (Note: only one header will be written; there will not be
        a header for each InstanceCatalog in the CompoundInstanceCatalog; default True)

        @param [in] write_mode is 'w' if you want to overwrite the output file or
        'a' if you want to append to an existing output file (default: 'w')
        """

        instantiated_ic_list = [None]*len(self._ic_list)

        # first, loop over all of the InstanceCatalog and CatalogDBObject classes, pre-processing
        # them (i.e. verifying that they have access to all of the columns they need)
        for ix, (icClass, dboClass) in enumerate(zip(self._ic_list, self._dbo_list)):
            best_connection = self.find_a_connection(dboClass)
            if best_connection is None:
                dbo = dboClass()
                self._active_connections.append(dbo.connection)
            else:
                dbo = dboClass(connection=best_connection)

            ic = icClass(dbo, obs_metadata=self._obs_metadata)
            ic._write_pre_process()
            instantiated_ic_list[ix] = ic

        for row in self._dbObjectGroupList:
            if len(row) == 1:
                ic = instantiated_ic_list[row[0]]
                ic._query_and_write(filename, chunk_size=chunk_size,
                                    write_header=write_header, write_mode=write_mode,
                                    obs_metadata=self._obs_metadata,
                                    constraint=self._constraint)
                write_mode = 'a'
                write_header = False

        default_compound_dbo = None
        if self._compoundDBclass is not None:
            if not hasattr(self._compoundDBclass, '__getitem__'):
                default_compound_dbo = CompoundCatalogDBObject
            else:
                for dbo in self._compoundDBclass:
                    if dbo._table_restriction is None:
                        default_compound_dbo = dbo
                        break

                if default_compound_dbo is None:
                    default_compound_dbo is CompoundCatalogDBObject

        for row in self._dbObjectGroupList:
            if len(row) > 1:
                dbObjClassList = [self._dbo_list[ix] for ix in row]
                catList = [instantiated_ic_list[ix] for ix in row]

                # if a connection is already open to the database, use
                # it rather than opening a new connection
                best_connection = self.find_a_connection(dbObjClassList[0])

                if self._compoundDBclass is None:
                    compound_dbo = CompoundCatalogDBObject(dbObjClassList, connection=best_connection)
                elif not hasattr(self._compoundDBclass, '__getitem__'):
                    # if self._compoundDBclass is not a list
                    try:
                        compound_dbo = self._compoundDBclass(dbObjClassList)
                    except:
                        compound_dbo = default_compound_dbo(dbObjClassList)
                else:
                    compound_dbo = None
                    for candidate in self._compoundDBclass:
                        use_it = True
                        if False in [candidate._table_restriction is not None
                                     and dbo.tableid in candidate._table_restriction
                                     for dbo in dbObjClassList]:

                            use_it = False

                        if use_it:
                            compound_dbo = candidate(dbObjClassList)
                            break

                    if compound_dbo is None:
                        compound_dbo = default_compound_dbo(dbObjClassList)

                self._write_compound(catList, compound_dbo, filename,
                                     chunk_size=chunk_size, write_header=write_header,
                                     write_mode=write_mode)
                write_mode = 'a'
                write_header = False

    def _write_compound(self, catList, compound_dbo, filename,
                        chunk_size=None, write_header=False, write_mode='a'):
        """
        Write out a set of InstanceCatalog instantiations that have been
        determined to query the same database table.

        @param [in] catList is the list of InstanceCatalog instantiations

        @param [in] compound_db is the CompoundCatalogDBObject instantiation
        associated with catList

        @param [in] filename is the name of the file to be written

        @param [in] chunk_size is an optional parameter telling the CompoundInstanceCatalog
        to query the database in manageable chunks (in case returning the whole catalog
        takes too much memory)

        @param [in] write_header a boolean specifying whether or not to add a header
        to the output catalog (Note: only one header will be written; there will not be
        a header for each InstanceCatalog in the CompoundInstanceCatalog; default True)

        @param [in] write_mode is 'w' if you want to overwrite the output file or
        'a' if you want to append to an existing output file (default: 'w')
        """

        colnames = []
        master_colnames = []
        name_map = []
        dbObjNameList = [db.objid for db in compound_dbo._dbObjectClassList]
        for name, cat in zip(dbObjNameList, catList):
            localNames = []
            local_map = {}
            for colName in cat._active_columns:
                colnames.append('%s_%s' % (name, colName))
                localNames.append('%s_%s' % (name, colName))
                local_map['%s_%s' % (name, colName)] = colName
            master_colnames.append(localNames)
            name_map.append(local_map)

        master_results = compound_dbo.query_columns(colnames=colnames,
                                                    obs_metadata=self._obs_metadata,
                                                    constraint=self._constraint,
                                                    chunk_size=chunk_size)

        with open(filename, write_mode) as file_handle:
            if write_header:
                catList[0].write_header(file_handle)

            new_dtype_list = [None]*len(catList)

            first_chunk = True
            for chunk in master_results:
                for ix, (catName, cat) in enumerate(zip(dbObjNameList, catList)):

                    if first_chunk:
                        for iy, name in enumerate(master_colnames[ix]):
                            if name not in chunk.dtype.fields:
                                master_colnames[ix][iy] = name_map[ix][name]

                    local_recarray = chunk[master_colnames[ix]].view(numpy.recarray)

                    local_recarray.flags['WRITEABLE'] = False  # so numpy does not raise a warning
                    # because it thinks we may accidentally
                    # write to this array
                    if new_dtype_list[ix] is None:
                        new_dtype = numpy.dtype([
                                                tuple([dd.replace(catName+'_', '')] +
                                                      [local_recarray.dtype[dd]])
                                                for dd in master_colnames[ix]
                                                ])
                        new_dtype_list[ix] = new_dtype

                    local_recarray.dtype = new_dtype_list[ix]
                    cat._write_recarray(local_recarray, file_handle)

                first_chunk = False
