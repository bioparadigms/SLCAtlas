#!/usr/bin/env python

import mysql.connector
import sql
import sql.aggregate
import logging
import collections.abc
import sqlite3

#logging.basicConfig(level=logging.DEBUG)

class Database(object):
    def __init__(self, host='127.0.0.1', database='your_db_name',
                 user='your_username', password='your_password'):
        self._name = database
        self._host = host
        self._user = user
        self._password = password
        self.connect()
        self.logger = logging.getLogger('mysql')
        #self.logger.setLevel(logging.DEBUG)

    def connect(self):
        self.conn = mysql.connector.connect(host=self._host,
                                            database=self._name,
                                            user=self._user,
                                            password=self._password)

    def commit(self):
        self.logger.debug('committing records')
        self.conn.commit()

    def close(self):
        self.conn.close()

    def query(self, *args, cursor=None, dictionary=False):
        try:
            if cursor is None: cur = self.conn.cursor(dictionary=dictionary)
            else: cur = cursor
            if isinstance(args[0], sql.Query):
                self.logger.debug(repr(tuple(args[0])))
                cur.execute(*tuple(args[0]))
            elif isinstance(args[0], str):
                self.logger.debug(repr(args))
                cur.execute(args[0], tuple(args[1:]))
            else:
                self.logger.debug(repr(args))
                cur.execute(*args)
        except (mysql.connector.errors.OperationalError, AttributeError):
            self.connect()
            if cursor is None: cur = self.conn.cursor(dictionary=dictionary)
            else: cur = cursor
            if isinstance(args[0], sql.Query):
                self.logger.debug(repr(tuple(args[0])))
                cur.execute(*tuple(args[0]))
            elif isinstance(args[0], str):
                self.logger.debug(repr(args))
                cur.execute(args[0], tuple(args[1:]))
            else:
                self.logger.debug(repr(args))
                cur.execute(*args)
        return cur

    def escape(s):
        # string escaping
        # https://dev.mysql.com/doc/refman/8.0/en/string-literals.html
        # % and _ need only be escaped in a LIKE context, " need not be escaped within '
        return s.translate(s.maketrans({
            '\0': '\\0',
            '\'': '\\\'',
            '\b': '\\b',
            '\n': '\\n',
            '\r': '\\r',
            '\t': '\\t',
            '\\': '\\\\',
            }))
    def format(v):
        if isinstance(v, int) or isinstance(v, bool): return '{:d}'.format(v)
        elif isinstance(v, float): return '{}'.format(v)
        elif v is None: return 'NULL'
        else: return '\'{}\''.format(Database.escape(str(v)))

class SQLiteDatabase(Database):
    def __init__(self, filename):
        super(SQLiteDatabase, self).__init__(database=filename, host=None, user=None, password=None)
        self.logger = logging.getLogger('sqlite')
        sql.Flavor.set(sql.Flavor(paramstyle='qmark'))
    def connect(self):
        self.conn = sqlite3.connect(self._name)
        self.conn.row_factory = sqlite3.Row
    def query(self, *args, cursor=None, dictionary=False):
        if cursor is None: cur = self.conn.cursor()
        else: cur = cursor
        if isinstance(args[0], sql.Query):
            self.logger.debug(repr(tuple(args[0])))
            cur.execute(*tuple(args[0]))
        elif isinstance(args[0], str):
            self.logger.debug(repr(args))
            cur.execute(args[0], tuple(args[1:]))
        else:
            self.logger.debug(repr(args))
            cur.execute(*args)
        return cur

from sql.functions import Function
class Substring_index(Function):
    __slots__ = ()
    _function = 'SUBSTRING_INDEX'

class Table(sql.Table):
    # "database" in the next line is not the same as schema
    def __init__(self, name, db, database=None, explicit_record_ids=True):
        super(Table, self).__init__(name, None, database)
        self._last_record_id = None
        # check if we have a record_id and a manual field
        info_schema = sql.Table('columns', 'information_schema')
        q = info_schema.select(
            info_schema.column_name,
            where=((info_schema.table_schema == db._name) &
                   (info_schema.table_name == name)))
        cur = db.query(q)
        self._has_record_id = False
        self._has_manual = False
        for row in cur.fetchall():
            if row[0] == 'record_id': self._has_record_id = True
            elif row[0] == 'manual': self._has_manual = True
        cur.close()
        self._explicit_record_ids = explicit_record_ids
        self._db = db

    def get_next_record_id(self):
        if not self._has_record_id: return None
        if self._last_record_id is not None:
            record_id = self._last_record_id
        else:
            q = self.select(sql.aggregate.Max(self.record_id))
            cur = self._db.query(q)
            record_id = cur.fetchone()[0]
            cur.close()
            if record_id is None: record_id = 0
        self._last_record_id = record_id + 1
        return record_id+1

    def print_insert_sql(self, row, manual=None, user=None, comments=None):
        keys = []
        values = []
        if isinstance(row, dict):
            for key in row.keys():
                keys.append('`{}`'.format(Database.escape(key)))
                values.append(Database.format(row[key]))
        elif isinstance(row, list) or isinstance(row, tuple):
            for key, value in row:
                keys.append('`{}`'.format(Database.escape(key)))
                values.append(Database.format(row[key]))
        if keys:
            if user is not None:
                keys.append('`user`')
                values.append(Database.format(user))
            if isinstance(comments, list) or isinstance(comments, tuple):
                keys.append('`comments`')
                values.append(Database.format('\n'.join(comments)))
            elif isinstance(comments, str):
                keys.append('`comments`')
                values.append(Database.format(comments))
            if self._has_manual and (manual is not None):
                keys.insert(0, '`manual`')
                values.insert(0, '{:d}'.format(manual))
            if self._has_record_id:
                keys.insert(0, '`record_id`')
                if self._explicit_record_ids:
                    values.insert(0, '{:d}'.format(self.get_next_record_id()))
                else:
                    values.insert(0, '@next_record_id')
                    print('SELECT MAX(`record_id`)+1 INTO @next_record_id '
                          'FROM {};'.format(Database.escape(str(self))))
            print('INSERT INTO {} ({}) VALUES ({});'.format(
                Database.escape(str(self)),
                ', '.join(keys),
                ', '.join(values)
            ))

class HistoryTable(sql.Table):
    def __init__(self, name, db, schema=None, database=None,
                 explicit_record_ids=True):
        super(HistoryTable, self).__init__(name, schema, database)
        self._history = sql.Table(name + '_history', schema, database)
        self._last_record_id = None
        # check if we have a manual field
        # record_id is mandatory here
        info_schema = sql.Table('columns', 'information_schema')
        q = info_schema.select(
            info_schema.column_name,
            where=((info_schema.table_schema == db._name) &
                   (info_schema.table_name == name)))
        cur = db.query(q)
        self._has_record_id = True
        self._has_manual = False
        for row in cur.fetchall():
            if row[0] == 'manual': self._has_manual = True
        cur.close()
        self._explicit_record_ids = explicit_record_ids
        self._db = db

    def get_next_record_id(self):
        if self._last_record_id is not None:
            record_id = self._last_record_id
        else:
            q = self._history.select(sql.aggregate.Max(self._history.record_id))
            cur = self._db.query(q)
            record_id = cur.fetchone()[0]
            cur.close()
            if record_id is None: record_id = 0
        self._last_record_id = record_id + 1
        return record_id+1

    def print_insert_sql(self, row, manual=0, user=None, comments=None):
        keys = []
        values = []
        if isinstance(row, dict):
            for key in row.keys():
                keys.append('`{}`'.format(Database.escape(key)))
                values.append(Database.format(row[key]))
        elif isinstance(row, list) or isinstance(row, tuple):
            for key, value in row:
                keys.append('`{}`'.format(Database.escape(key)))
                values.append(Database.format(row[key]))
        if keys:
            if user is not None:
                keys.append('`user`')
                values.append(Database.format(user))
            if isinstance(comments, list) or isinstance(comments, tuple):
                keys.append('`comments`')
                values.append(Database.format('\n'.join(comments)))
            elif isinstance(comments, str):
                keys.append('`comments`')
                values.append(Database.format(comments))
            if self._has_manual:
                keys.insert(0, '`manual`')
                values.insert(0, '{:d}'.format(manual))
            keys.insert(0, '`record_id`')
            if self._explicit_record_ids:
                values.insert(0, '{:d}'.format(self.get_next_record_id()))
            else:
                values.insert(0, '@next_record_id')
                print('SELECT MAX(`record_id`)+1 INTO @next_record_id '
                      'FROM {};'.format(Database.escape(str(self._history))))
            print('INSERT INTO {} ({}) VALUES ({});'.format(
                Database.escape(str(self._history)),
                ', '.join(keys),
                ', '.join(values)
            ))

    def print_update_sql(self, _id, values, manual=0, user=None, comments=None):
        # write sql query
        print('CREATE TEMPORARY TABLE `tmptable` '
              'SELECT * FROM {} WHERE `id` = {:d}; '
              'UPDATE `tmptable` '
              'SET `id` = NULL;'.format(str(self._history), _id))
        # modify
        if isinstance(comments, list) or isinstance(comments, tuple):
            c = Database.format('\n'.join(comments))
        elif isinstance(comments, str):
            c = Database.format(comments)
        print('UPDATE `tmptable` '
              'SET `user` = {}, '
              '`comments` = {};'.format(Database.format(user), c))
        if self._has_manual:
            print('UPDATE `tmptable` SET `manual` = {:d};'.format(manual))
        for key, value in values.items():
            print('UPDATE `tmptable` SET `{}` = {};'.format(
                Database.escape(key), Database.format(value)))
        print('INSERT INTO {} SELECT * FROM `tmptable`; '
              'DROP TEMPORARY TABLE IF EXISTS '
              '`tmptable`;'.format(str(self._history)))

    """
        db: Database() instance
        history_table: sql.Table() instance
        field: string
        key: dict -> {field: value} will be used to index history_table
    """
    def get_field_changes(self, fields, key):
        cond = lambda table: ' AND '.join(
            ['{!s}.`{!s}` = {}'.format(table, k, Database.format(v))
             for k, v in key.items()])
        select_fields = []
        changed_conds = []
        if isinstance(fields, str):
            select_fields.append("t2.`{!s}` AS old_value".format(fields))
            select_fields.append("t1.`{!s}` AS new_value".format(fields))
            changed_conds.append("(BINARY t1.`{0!s}` != t2.`{0!s}`)".format(fields))
        elif isinstance(fields, collections.abc.Iterable):
            for i, field in enumerate(fields):
                select_fields.append("t2.`{!s}` AS old_value{:d}".format(field, i+1))
                select_fields.append("t1.`{!s}` AS new_value{:d}".format(field, i+1))
                changed_conds.append("(BINARY t1.`{0!s}` != t2.`{0!s}`)".format(field))
        else:
            TypeError('"fields" needs to be string or iterable')
        changed_conds.append('(t1.`deleted` != t2.`deleted`)')
        ## TODO: make the "manual" field optional based on self._has_manual
        q = """\
            SELECT
                t1.id,
                t1.record_id,
                {fields},
                t1.manual,
                t1.deleted,
                t1.timestamp,
                IFNULL({conditions}, 1) AS changed
            FROM {table} AS t1
            LEFT JOIN {table} AS t2
            ON t2.id = (
                SELECT
                    MAX(t3.id)
                FROM {table} AS t3
                WHERE ({cond_t3}) AND (t3.id < t1.id)
                )
            WHERE {cond_t1}
            HAVING changed = 1
            ORDER BY t1.id DESC;""".format(
                table=str(self._history),
                fields=', '.join(select_fields),
                conditions=' OR '.join(changed_conds),
                cond_t1=cond('t1'),
                cond_t3=cond('t3'))
        cur = self._db.query(q, dictionary=True)

        #cur = self._db.query("""\
        #               SELECT
        #                   t1.id,
        #                   t1.record_id,
        #                   t2.{field} AS old_value,
        #                   t1.{field} AS new_value,
        #                   t1.manual,
        #                   t1.deleted,
        #                   IFNULL((BINARY t1.{field} != t2.{field}) OR 
        #                       (t1.deleted != t2.deleted), 1) AS changed
        #               FROM {table} AS t1
        #               LEFT JOIN {table} AS t2
        #               ON t2.id = (
        #                   SELECT
        #                       MAX(t3.id)
        #                   FROM {table} AS t3
        #                   WHERE ({cond_t3}) and (t3.id < t1.id)
        #                   )
        #               WHERE {cond_t1}
        #               HAVING changed = 1
        #               ORDER BY t1.id DESC;
        #               """.format(table=str(self._history),
        #                          field='`{}`'.format(field),
        #                          cond_t3=cond('t3'),
        #                          cond_t1=cond('t1')),
        #               dictionary=True)
        changes = cur.fetchall()
        cur.close()
        if not isinstance(fields, str) and isinstance(fields, collections.abc.Iterable):
            for change in changes:
                change['old_value'] = tuple(change['old_value{:d}'.format(x)]
                                            for x in range(1, len(fields)+1))
                change['new_value'] = tuple(change['new_value{:d}'.format(x)]
                                            for x in range(1, len(fields)+1))
                for x in range(1, len(fields)+1):
                    del change['old_value{:d}'.format(x)]
                    del change['new_value{:d}'.format(x)]
        return changes

    def insert_row(self, new_row, manual=0, keys=None, user=None,
                   comments=None):
        if keys is None: keys = new_row.keys()
        row = dict()
        ## TODO: make "manual" field optional based on self._has_manual
        #row.append(('manual', manual))
        #row['manual'] = manual
        for key in keys:
            if key in new_row:
                row[key] = new_row[key]
        self.print_insert_sql(row, manual=manual, user=user, comments=comments)

    # if we are updating a field that has been changed manually, don't update if
    # the update would revert to the state before the manual change.
    # otherwise, update even if the field has been manually changed.
    def update_row(self, old_row, to_update, manual=0, keys=None, user=None,
                   comments=None):
        changed_fields = []
        if keys is None: keys = to_update.keys()
        for key in keys:
            if (key not in old_row) or (to_update[key] != old_row[key]):
                changed_fields.append(key)
        ## TODO: make "manual" field optional everywhere here based in
        ## self._has_manual
        # if there are changes, check the history if any of those fields have been
        # set manually
        # POLICY HERE: we don't revert to the last automatic value if the most
        # recent changes were manual.
        for field in changed_fields[:]:
            # get a list of change events
            changes = self.get_field_changes(field, {'record_id':
                                                     old_row['record_id']})

            # look for last automatic change
            for change in changes:
                if change['manual'] == 0: break
            else:
                # no automatic changes found?
                # then, don't modify!
                changed_fields.remove(field)
                continue
            # if last change in history is automatic, then
            # change['new_value'] == old_row[field] != to_update[field]
            # so the next if will be false
            # and field will be changed -> OK
            if change['new_value'] == to_update[field]:
                # don't revert
                logging.warning('not reverting to manually overridden value '
                                '{!r} for field {} in record # {:d} in table '
                                '{!s}'.format(to_update[field], field,
                                              old_row['record_id'],
                                              self._history))

                changed_fields.remove(field)

        if changed_fields:
            self.print_update_sql(old_row['id'],
                                  dict([(k, v) 
                                        for k, v in to_update.items() 
                                        if k in changed_fields]),
                                  manual=manual, user=user, comments=comments)

    """
        db: Database() instance
        history_table: sql.Table() instance
        fields: string or list (iterable)
        key: dict -> see get_field_changes
        new_values: dict (new values, ALL of them)
    """
    def update_one_to_many(self, fields, key, new_values, manual=0,
                           user=None, comments=None):
        #cond = lambda table: ' AND '.join(
        #    ['{!s}.`{}` = {}'.format(table, k, Database.format(v))
        #     for k, v in key.items()])
        q = ('SELECT DISTINCT `record_id` FROM {!s} '
             'WHERE {}'.format(
                 self._history, ' AND '.join(
                     ['{!s}.`{}` = {}'.format(str(self._history), k, Database.format(v))
                      for k, v in key.items()])))
        cur = self._db.query(q)
        record_ids = [x[0] for x in cur.fetchall()]
        cur.close()

        # we collect records here that can be re-used for insertion
        reusables = dict([(x, []) for x in new_values])
        to_delete = []
        # get change events
        # decide fate of each record
        for record_id in record_ids:
            changes = self.get_field_changes(fields,
                                             {'record_id': record_id})
            if changes[0]['new_value'] in new_values:
                if changes[0]['deleted'] == 0:
                    # the value is already live
                    new_values.remove(changes[0]['new_value'])
                    logging.debug('removing {!r} from '
                                  'new_values, new_values after = '
                                  '{!r}'.format(changes[0]['new_value'],
                                                new_values))
                else:
                    if changes[0]['manual'] == 1:
                        # the value was manually deleted, don't revert
                        # unless this is a manual update!
                        if manual == 0:
                            new_values.remove(changes[0]['new_value'])
                            logging.warning(
                                'value {} for key {!r} field {!r} in table '
                                '{!s} was manually deleted, not '
                                'reverting.'.format(
                                    changes[0]['new_value'], key, fields,
                                    self._history))
                        else:
                            logging.warning(
                                'undeleting manually deleted value {} '
                                'for key {!r} field {!r} in table {!s}'.format(
                                    changes[0]['new_value'], key, fields,
                                    self._history))
                            self.print_update_sql(
                                changes[0]['id'], {'deleted': 0},
                                manual=manual, user=user, comments=comments)
                    else:
                        # the value was automatically deleted, undelete
                        logging.warning(
                            'undeleting automatically deleted value {} '
                            'for key {!r} field {!r} in table {!s}'.format(
                                changes[0]['new_value'], key, fields,
                                self._history))
                        self.print_update_sql(
                            changes[0]['id'], {'deleted': 0},
                            manual=manual, user=user, comments=comments)
                        new_values.remove(changes[0]['new_value'])
            else:
                # find last automatic value
                for change in changes:
                    if change['manual'] == 0: break
                if ((change['manual'] == 0) and 
                    (change['new_value'] in new_values)):
                    # value is still being overridden, don't revert
                    # unless this is a manual update!
                    if manual == 0:
                        logging.warning(
                            'not reverting to overridden value {} for '
                            'key {!r} field {!r} in table {!s}'.format(
                                change['new_value'], key, fields,
                                self._history))
                        new_values.remove(change['new_value'])
                    else:
                        logging.warning(
                            'reverting to overridden value {} for '
                            'key {!r} field {!r} in table {!s}'.format(
                                change['new_value'], key, fields,
                                self._history))
                        if isinstance(fields, str):
                            self.print_update_sql(
                                changes[0]['id'], {field: change['new_value']},
                                manual=manual, user=user, comments=comments)
                        elif isinstance(fields, collections.abc.Iterable):
                            self.print_update_sql(
                                changes[0]['id'], 
                                dict(zip(fields, change['new_value'])),
                                manual=manual, user=user, comments=comments)
                        else:
                            TypeError('"fields" must be string or '
                                      'iterable')
                        new_values.remove(change['new_value'])
                else:
                    # we could delete this record if still alive,
                    # or we could re-use it.
                    # also, don't delete records that were initially manual
                    # unless this is a manual update!
                    if ((changes[0]['deleted'] == 0) and
                        (change['manual'] == 0)):
                        to_delete.append((change['record_id'],
                                          changes[0]['id']))
                    if manual == 0:
                        if ((changes[0]['deleted'] == 0) and
                            (change['manual'] == 1)):
                            logging.warning(
                                'not deleting manually set value {} for '
                                'key {!r} field {!r} in table {!s}'.format(
                                    changes[0]['new_value'], key, fields,
                                    self._history))
                    else:
                        logging.warning(
                            'deleting record #{:d} manually with value '
                            '{} for key {!r} field {!r} in table '
                            '{!s}'.format(
                                change['record_id'],
                                changes[0]['new_value'], key, fields,
                                self._history))
                        to_delete.append((change['record_id'],
                                          changes[0]['id']))
                    # here, we check the whole history whether it can be
                    # re-used
                    for change in changes:
                        if change['new_value'] in new_values:
                            reusables[change['new_value']].append(
                                (change['record_id'], changes[0]['id']))
        # these are values that we need to insert
        for value in new_values:
            if reusables[value]:
                # there are some records we could re-use
                record_id, _id = reusables[value].pop(0)
                # re-use record, undelete if necessary
                logging.warning('re-using record #{:d} for key {!r} field '
                                '{!r} in table {!s}'.format(
                                    record_id, key, fields, self._history))
                if isinstance(fields, str):
                    self.print_update_sql(_id, {fields: value, 'deleted': 0},
                                          manual=manual, user=user,
                                          comments=comments)
                elif isinstance(fields, collections.abc.Iterable):
                    self.print_update_sql(_id, {**dict(zip(fields, value)),
                                                **{'deleted': 0}},
                                          manual=manual, user=user,
                                          comments=comments)

                # don't delete this record since we re-used it
                to_delete.remove(record_id)
            else:
                # insert new record
                #print_insert_sql(...)
                logging.info('new value {} for key {!r} field {!r} in table '
                             '{!s}'.format(value, key, fields, self._history))
                if isinstance(fields, str):
                    self.print_insert_sql({**key, **{fields: value}},
                                          manual=manual,
                                          user=user, comments=comments)
                elif isinstance(fields, collections.abc.Iterable):
                    self.print_insert_sql({**key,
                                           **dict(zip(fields, value))},
                                          manual=manual,
                                          user=user, comments=comments)
        # then delete all remaining records
        for record_id, _id in to_delete:
            logging.info('deleting record #{:d} for key {!r} field {!r} in '
                         'table {!s}'.format(
                             record_id, key, fields, self._history))
            self.print_update_sql(_id, {'deleted': 1}, manual=manual,
                                  user=user, comments=comments)

