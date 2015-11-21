import csv
import os

# Import psycopg2 if possible, otherwise data access will only be available from local files
try:
    import psycopg2
except:
    print('psycopg2 is unavailable, no database access will be possible')
    
# Initialize the database connection, if possible
pg_database = os.environ.get('POSTGRES_CONNECTION_DATABASE')
pg_user = os.environ.get('POSTGRES_CONNECTION_USER')
if pg_database is not None and pg_user is not None:
  conn = psycopg2.connect(database=pg_database, user=pg_user)
  conn.autocommit = True

# Pulls all columns from the specified database table or local csv file
def get_bounded_data(table_name, bounds):
  if 'conn' in globals():
      cursor = conn.cursor()
      cursor.execute('SELECT * FROM %s WHERE (lat BETWEEN %s AND %s) and (lon BETWEEN %s AND %s)' % (table_name, bounds['lat_min'], bounds['lat_max'], bounds['lon_min'], bounds['lon_max']))
      return cursor.fetchall()
  else:
      reader = csv.reader(open('test_data/%s.csv' % table_name, 'rb'))
      return list(reader)
