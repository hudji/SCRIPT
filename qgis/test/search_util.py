from __future__ import print_function, division, absolute_import

import json
import time
import requests
import landsat
import landsat.utils
from landsat import *
from landsat.utils import three_digit, create_paired_list, geocode
from landsat import settings
import geocoder

class Search(object):
    """ The search class """
    def __init__(self):
        self.api_url = settings.API_URL
        def search(self, paths_rows=None, lat=None, lon=None, address=None, start_date=None, end_date=None,
                   cloud_min=None, cloud_max=None, limit=1, geojson=False):
            search_string = self.query_builder(paths_rows, lat, lon, address, start_date, end_date, cloud_min, cloud_max)

            # Have to manually build the URI to bypass requests URI encoding
            # The api server doesn't accept encoded URIs

            r = requests.get('%s?search=%s&limit=%s' % (self.api_url, search_string, limit))
            r_dict = json.loads(r.text)
            result = {}

            if 'error' in r_dict:
                result['status'] = u'error'
                result['code'] = r_dict['error']['code']
                result['message'] = r_dict['error']['message']

            elif 'meta' in r_dict:
                if geojson:
                    result = {
                        'type': 'FeatureCollection',
                        'features': []
                    }
                    for r in r_dict['results']:
                        feature = {
                            'type': 'Feature',
                            'properties': {
                                'sceneID': r['sceneID'],
                                'row': three_digit(r['row']),
                                'path': three_digit(r['path']),
                                'thumbnail': r['browseURL'],
                                'date': r['acquisitionDate'],
                                'cloud': r['cloudCoverFull']
                            },
                            'geometry': {
                                'type': 'Polygon',
                                'coordinates': [
                                    [
                                        [r['upperLeftCornerLongitude'], r['upperLeftCornerLatitude']],
                                        [r['lowerLeftCornerLongitude'], r['lowerLeftCornerLatitude']],
                                        [r['lowerRightCornerLongitude'], r['lowerRightCornerLatitude']],
                                        [r['upperRightCornerLongitude'], r['upperRightCornerLatitude']],
                                        [r['upperLeftCornerLongitude'], r['upperLeftCornerLatitude']]
                                    ]
                                ]
                            }
                        }

                        result['features'].append(feature)
                    else:
                        result['status'] = u'SUCCESS'
                        result['total'] = r_dict['meta']['results']['total']
                        result['limit'] = r_dict['meta']['results']['limit']
                        result['total_returned'] = len(r_dict['results'])
                        result['results'] = [{'sceneID': i['sceneID'],
                                              'sat_type': u'L8',
                                              'path': three_digit(i['path']),
                                              'row': three_digit(i['row']),
                                              'thumbnail': i['browseURL'],
                                              'date': i['acquisitionDate'],
                                              'cloud': i['cloudCoverFull']}
                                             for i in r_dict['results']]

                    return result

                    def query_builder(self, paths_rows=None, lat=None, lon=None, address=None, start_date=None, end_date=None,cloud_min=None, cloud_max=None):
                        query = []
                        or_string = ''
                        and_string = ''
                        search_string = ''

                        if paths_rows:
                            # Coverting rows and paths to paired list
                            new_array = landsat.create_paired_list(paths_rows)
                            paths_rows = ['(%s)' % self.row_path_builder(i[0], i[1]) for i in new_array]
                            or_string = '+OR+'.join(map(str, paths_rows))

                        if start_date and end_date:
                            query.append(self.date_range_builder(start_date, end_date))
                        elif start_date:
                            query.append(self.date_range_builder(start_date, '2100-01-01'))
                        elif end_date:
                            query.append(self.date_range_builder('2009-01-01', end_date))

                        if cloud_min and cloud_max:
                            query.append(self.cloud_cover_prct_range_builder(cloud_min, cloud_max))
                        elif cloud_min:
                            query.append(self.cloud_cover_prct_range_builder(cloud_min, '100'))
                        elif cloud_max:
                            query.append(self.cloud_cover_prct_range_builder('-1', cloud_max))

                        if address:
                            query.append(self.address_builder(address))
                        elif (lat is not None) and (lon is not None):
                            query.append(self.lat_lon_builder(lat, lon))

                        if query:
                            and_string = '+AND+'.join(map(str, query))

                        if and_string and or_string:
                            search_string = and_string + '+AND+(' + or_string + ')'
                        else:
                            search_string = or_string + and_string

                        return search_string

                    def row_path_builder(self, path='', row=''):

                        return 'path:%s+AND+row:%s' % (path, row)

                        def date_range_builder(self, start='2013-02-11', end=None):

                            if not end:
                                end = time.strftime('%Y-%m-%d')

                            return 'acquisitionDate:[%s+TO+%s]' % (start, end)

                    def cloud_cover_prct_range_builder(self, min=0, max=100):

                        return 'cloudCoverFull:[%s+TO+%s]' % (min, max)

                    def address_builder(self, address):

                        geocoded = landsat.geocode(address)
                        return self.lat_lon_builder(**geocoded)

                    def lat_lon_builder(self, lat=0, lon=0):

                        return ('upperLeftCornerLatitude:[%s+TO+1000]+AND+lowerRightCornerLatitude:[-1000+TO+%s]'
                                '+AND+lowerLeftCornerLongitude:[-1000+TO+%s]+AND+upperRightCornerLongitude:[%s+TO+1000]'
                                % (lat, lat, lon, lon))

#search = Search()
#search('003,003', '2014-01-01', '2014-06-01')
search = Search('003,003', '2014-01-01', '2014-06-01')