"""
Script to correct the WCS information from fits images using astrometry.net solving 
"""

from __future__ import print_function

import re
import os
import sys
import urllib as url
from urllib2 import urlopen
from urllib2 import Request
from urllib2 import HTTPError
from urllib import urlencode
from urllib import quote
import time

from client import Client

DIR = #Directory to load the images

def read_files():
    with open(DIR+'fitslist', 'r') as fp:
        files = fp.read().splitlines()
    return files

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--server', dest='server', default='http://nova.astrometry.net/api/',
                      help='Set server base URL (eg, http://nova.astrometry.net/api/)')
    parser.add_option('--apikey', '-k', dest='apikey',
                      help='API key for Astrometry.net web service; if not given will check AN_API_KEY environment variable')
    parser.add_option('--newfits', dest='newfits', help='Download resulting new-image.fits file, saving to given filename; implies --wait if --urlupload or --upload')
    parser.add_option('--wait', '-w', dest='wait', action='store_true', help='After submitting, monitor job status')
    parser.add_option('--jobid', '-i', dest='solved_id', type=int,help='retrieve result for jobId instead of submitting new image')
    parser.add_option('--substatus', '-s', dest='sub_id', help='Get status of a submission')
    parser.add_option('--wcs', dest='wcs', help='Download resulting wcs.fits file, saving to given filename; implies --wait if --urlupload or --upload')

    #parser.add_option('--newfits', dest='newfits', help='Download resulting new-image.fits file, saving to given filename; implies --wait if --urlupload or --upload')

    opt,args = parser.parse_args()
    if opt.apikey is None:
        # try the environment
        opt.apikey = os.environ.get('AN_API_KEY', None)
    if opt.apikey is None:
        parser.print_help()
        print()
        print('You must either specify --apikey or set AN_API_KEY')
        sys.exit(-1)
    useclient = True
    if useclient:
        client = Client(apiurl=opt.server)
        client.login(opt.apikey)

    files = read_files() 

    for filenames in files:
        opt.solved_id = None
        filename = DIR+filenames
        if useclient:
            upres = client.upload(filename)
            print(client.submission_images(1))
            opt.sub_id = upres['subid']
        else:
            cmd = "python client.py --server %s --apikey %s --upload \"%s\"" % (opt.server, opt.apikey, filename)
            print(cmd)
            os.system(cmd)

        if opt.wait:
            if opt.solved_id is None:
                if opt.sub_id is None:
                    print("Can't --wait without a submission id or job id!")
                    sys.exit(-1)

                while True:
                    stat = client.sub_status(opt.sub_id, justdict=True)
                    print('Got status:', stat)
                    jobs = stat.get('jobs', [])
                    if len(jobs):
                        for j in jobs:
                            if j is not None:
                                break
                        if j is not None:
                            print('Selecting job id', j)
                            opt.solved_id = j
                            break
                    time.sleep(5)

            while True:
                stat = client.job_status(opt.solved_id, justdict=True)
                print('Got job status:', stat)
                if stat.get('status','') in ['success']:
                    success = (stat['status'] == 'success')
                    break
                time.sleep(5)

        if opt.solved_id:
        # we have a jobId for retrieving results
            retrieveurls = []
            if opt.wcs:
                # We don't need the API for this, just construct URL
                url = opt.server.replace('/api/', '/wcs_file/%i' % opt.solved_id)
                retrieveurls.append((url, opt.wcs))
            if opt.newfits:
                url = opt.server.replace('/api/', '/new_fits_file/%i/' % opt.solved_id)
                retrieveurls.append((url, opt.newfits))

            for url,fn in retrieveurls:
                print('Retrieving file from', url, 'to', fn)
                f = urlopen(url)
                txt = f.read()
                fn = DIR+'WCS_CORRECTED/'+filenames
                w = open(fn, 'wb')
                w.write(txt)
                w.close()
                print('Wrote to', fn)
