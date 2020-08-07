#!/usr/bin/env python
""" STEP AddToDatabase- Version 1.0.0

    This step registers an image with the local image database.

    It does this by parsing the database configuration file, which's path is given in the 
    step config file, for the desired fields from FITS HDUs. It then gets field values for those 
    fields from self.datain and executes an SQL statement adding the record. 

    Assumptions necessary for the step to run without error:
    -the 'SEO' database and 'fits_data' table within it must both have been created 
    -the database config file, the path of which is in the step config, must exist and not be empty
    -the db config file's first entry, and thus the primary key for the database is 'file_path'
    -all database fields in the config file aside from 'file_path' map to FITS HDU fields 
    by the sql_field_to_hdu_field function below (will be true if they follow the convention in the config)
    
    @author: Enrique Collin
"""

import logging 
from darepype.drp import StepParent
import mysql.connector as mysql
from os import path

class StepAddToDatabase(StepParent):
    """ HAWC Pipeline Step Parent Object
        The object is callable. It requires a valid configuration input
        (file or object) when it runs.
    """
    stepver = '0.1' # pipe step version

    def parse_config(self, path_to_config):
        """
        Parses given database config file the SQL table fields within it. 
        
        Note that fields  are from the config file, not the database itself, 
        so if it has changed and the database has not been updated there will be an inconsistency there.

        Parameters:
        path_to_config: Config file for fits_data table with mysql database. Every line must be blank,
        start with #, or be of the format "<field name> <whitespace> <field_type (optional)>
        Note that this is the config for the database table, not for any pipestep
        
        Returns:
        Returns a list of SQL field names specified in the config file.
        Throws an error if the config file does not exist
        """
        sql_fields = []
        if not path.exists(path_to_config):
            err_msg = 'Path to database config given in StepAddToDatabase config does not exist!'
            self.log.error(err_msg)
            raise FileNotFoundError(err_msg)
        with open(path_to_config, 'r') as config:
            for line in config:
                trimmed = line.strip()
                # If line is a comment or only whitespace keep going
                if trimmed.startswith("#") or len(trimmed) == 0:
                    continue
                # Take only the first word of each line which is the field
                field = trimmed.split()[0]
                sql_fields.append(field)
        return sql_fields
    
    def sql_field_to_hdu_field(self, sql_field):
        """
        Given an sql_field that corresponds to an hdu field and follows
        naming convention, will convert it to the corresponding hdu_field.
        
        sql field naming convention mentioned above: take the name of the corresponding HDU field,
        make it lowercase, replace hyphens with underscores, and optionally append a trailing
        underscore (which is required if the application of the other rules turns the 
        HDU field into an an SQL keyword, like with "DEC"). 

        Parameters:
        sql_field: the string name of an sql_field that corresponds to a FITS HDU field
        The field name must follow naming convention which is: take the name of the HDU field,
        make it lowercase, replace hyphens with underscores, and

        Returns:
        Takes an sql field under the name convention described above and returns its HDU field
        counterpart.  
        """
        hdu_field = sql_field.upper().replace('_', '-')
        if hdu_field[-1] == '-':
            return hdu_field[:-1]
        else:
            return hdu_field
    
    def setup(self):
        """ ### Names and Prameters need to be Set Here ###
            Sets the internal names for the function and for saved files.
            Defines the input parameters for the current pipe step.
            Setup() is called at the end of __init__
            The parameters are stored in a list containing the following
            information:
            - name: The name for the parameter. This name is used when
                    calling the pipe step from command line or python shell.
                    It is also used to identify the parameter in the pipeline
                    configuration file.
            - default: A default value for the parameter. If nothing, set
                       '' for strings, 0 for integers and 0.0 for floats
            - help: A short description of the parameter.
        """
        ### Set Names
        # Name of the pipeline reduction step
        self.name='AddToDatabase'
        # Shortcut for pipeline reduction step and identifier for
        # saved file names.
        self.procname = 'unk'
        # Set Logger for this pipe step
        self.log = logging.getLogger('hawc.pipe.step.%s' % self.name)
        ### Set Parameter list
        # Clear Parameter list
        self.paramlist = []
        # Append parameters
        self.paramlist.append(['sql_username', 'default_user',
            'Username for the SQL user to access the database; only add priveleges needed'])
        self.paramlist.append(['sql_password', '',
            'Password for the SQL user to access the database'])
        self.paramlist.append(['database_config_path', './pipeconf_dbasereg.txt',
            'Path to database configuration file'])
        self.paramlist.append(['overwrite', False,
            ('If True, when StepAddToDatabase is run on a file that is already in '
            'the database, then it will be overwritten. Requires an SQL user with ' 
            'delete permissions. If overwrite is False and a duplicate file exists, '
            'an error will be thrown.')])

    def run(self):
        """ 
        Connects to the SEO database and adds self.datain to the database.
        Only adds fields described in the database config file, which's path is in the step
        config file. 

        Assumptions:
        -SEO database has been created
        -file_path is the first entry in the DB config file and thus primary key for the DB
        -The file being added is not already in the database. 
        If any of these is violated an error is thrown. 
        """
        # Copy datain to dataout (the data doesn't actually have to change)
        self.dataout = self.datain
        # The user should be a mySQL user granted ONLY add permissions
        SQL_user = self.getarg('sql_username')
        SQL_pass = self.getarg('sql_password')
        database_config_path = self.getarg('database_config_path')
        overwrite = self.getarg('overwrite')
        
        try:
            db = mysql.connect(
                host="localhost",
                user=SQL_user,
                passwd=SQL_pass,
                auth_plugin='mysql_native_password'
            )
        except mysql.errors.ProgrammingError as err:
            err_msg = (f'Encountered error connecting to DB as "{SQL_user}". '
                        'User or pass may be wrong')
            self.log.error(err_msg)
            raise RuntimeError(err_msg) from err
        self.log.info(f'Successfully connected to SQL server as "{SQL_user}"')        
        cursor = db.cursor()

        try:
            cursor.execute('USE seo;')
        except mysql.errors.ProgrammingError as err:
            err_msg = ('Encountered error running SQL command "USE seo". Could mean' 
                      'seo database does not exist for some reason. Perhaps '
                      'create_database.py needs to be run')
            self.log.error(err_msg)    
            cursor.close()
            db.close()
            raise RuntimeError(err_msg) from err

        sql_fields = self.parse_config(database_config_path)
        self.log.debug(f'Successfully read in database config at {database_config_path}')        
        if sql_fields[0] != 'file_path':
            err_msg = ('Current StepAddToDatabase code assumes the first entry'
                        'in the database config file is "file_path", but this '
                        'is not the case!')
            self.log.error(err_msg)    
            raise RuntimeError(err_msg)
        

        fields_str = ', '.join(sql_fields)
        # Get format specifier for each sql_field; truncate tailing space and comma
        format_specifiers = ('%s, ' * len(sql_fields))[:-2]
        insert_query = f'INSERT INTO fits_data ({fields_str}) VALUES ({format_specifiers});'

        datain_field_vals = []
        if not path.exists(self.datain.filename):
            self.log.warning(
                (f'Starting to add record for {self.datain.filename} to database'
                ' even though os.path.exists for it is False; it has not yet been saved')
            )

        datain_field_vals.append(self.datain.filename)
        for sql_field in sql_fields:
            if sql_field != 'file_path':
                hdu_field = self.sql_field_to_hdu_field(sql_field)
                # Need val in string form for the below.
                val = self.datain.header[hdu_field]
                if isinstance(val, bool):
                    datain_field_vals.append(int(val))
                else:
                    datain_field_vals.append(val)

        self.log.debug(
            ('About to attempt to execute the following SQL: '
            f'"{insert_query}" with values "{datain_field_vals}"')
        )

        if overwrite:
            self.log.debug('Duplicate encountered; about to attempt to delete it')
            
            delete_cmd = f'DELETE FROM fits_data WHERE file_path = "{self.datain.filename}";'
            try: 
                cursor.execute(delete_cmd)
            except mysql.errors.ProgrammingError as err:
                err_msg = 'Error trying to overwrite duplicate file; likely insufficient privileges'
                self.log.error(err_msg)
                cursor.close()
                db.close()
                raise RuntimeError(err_msg) from err
            self.log.info('Successfully deleted duplicate file from database')
   
        try:
            cursor.execute(insert_query, tuple(datain_field_vals))
            db.commit()
        except mysql.errors.ProgrammingError as err:
            err_msg = 'The error could mean the config file is not up to date with db'
            self.log.error(err_msg)
            cursor.close()
            db.close()
            raise RuntimeError(err_msg) from err
        except mysql.errors.IntegrityError as err:
            err_msg = "The error likely means the file you're adding to the db is already there"
            self.log.error(err_msg)
            cursor.close()
            db.close()
            raise RuntimeError(err_msg) from err
        
        self.log.info('Successfully added file %s to the database' % self.datain.filename)
        cursor.close()
        db.close()


    def test(self):
        """ Test Pipe Step Parent Object:
            Runs a adds to the database, checking for errors.
            Should not be run on the actual database as it will add bogus records
        """
        # log message
        self.log.info('Testing pipe step %s' %self.name)
        # log message
        self.log.info('Testing pipe step %s - Done' %self.name)
    
if __name__ == '__main__':
    """ Main function to run the pipe step from command line on a file.
        Command:
          python stepparent.py input.fits -arg1 -arg2 . . .
        Standard arguments:
          --config=ConfigFilePathName.txt : name of the configuration file
          -t, --test : runs the functionality test i.e. pipestep.test()
          --loglevel=LEVEL : configures the logging output for a particular level
          -h, --help : Returns a list of 
    """
    StepAddToDatabase().execute()


""" === History ===
August 3, 2020-First working version created. Assumes file_path is primary key and only non-HDU field.
"""
