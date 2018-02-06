#!/usr/bin/env python
#    This file is part of Acca plugin.

#    Acca plugin is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    Acca plugin is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Acca plugin.  If not, see <http://www.gnu.org/licenses/>.


from toar import *
from acca import *
from plugin import *
import optparse
import os

def main(metafile,tmpdir,maskfile,debuglevel=0,shadows=True,cloud_sign=True,single=False):
    t=CToar(metafile,tmpdir,debuglevel)
    t.start()
    t.wait()
    a=CAcca(os.path.join(tmpdir,os.path.basename(metafile)),maskfile,debuglevel,shadows,cloud_sign,single)
    a.start()
    a.wait()

if __name__=="__main__":
    parser=optparse.OptionParser()
    parser.add_option("--verbose","-v",default=False,action="store_true",help="Verbose output")
    parser.add_option("--shadows","-s",default=False,action="store_true",help="Shadows")
    parser.add_option("--cloudsign","-c",default=False,action="store_true",help="Cloud signature")
    parser.add_option("--singlepass","-p",default=False,action="store_true",help="Single path")
    parser.add_option("--metafile",dest="metafile",help="Path to metafile")
    parser.add_option("--maskfile",dest="maskfile",help="Path to maskfile")
    parser.add_option("--tmpdir",dest="tmpdir",help="Path to tmpdir")
    (options,args)=parser.parse_args()
    main(options.metafile,
        options.tmpdir,
        options.maskfile,
        1 if options.verbose else 0,
        options.shadows,
        options.cloudsign,
        options.singlepass)
