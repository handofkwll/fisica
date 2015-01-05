from __future__ import absolute_import

import collections
import datetime
import os
import shutil

from extern import mako
from mako.template import Template
from mako.lookup import TemplateLookup

import templates
import templates.resources as resources


class Renderer(object):
    """Class to render results into html.
    """

    def __init__(self, result):
        self.result = result

    def copy_resources(self, dirname):
        outdir = os.path.join(dirname, 'resources')

        # shutil.copytree complains if the output directory exists
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        # copy all uncompressed non-python resources to output directory
        src = os.path.dirname(resources.__file__)
        dst = outdir
        ignore_fn = shutil.ignore_patterns('*.zip','*.py','*.pyc','CVS*')
        shutil.copytree(src,
                        dst,
                        symlinks=False,
                        ignore=ignore_fn)

#        # unzip fancybox to output directory
#        infile = os.path.join(src, 'fancybox.zip')
#        z = zipfile.ZipFile(infile, 'r')
#        z.extractall(outdir)

    def run(self):
        # make a directory to contain this result
        now = datetime.datetime.today()        
        dirname = 'fisica-%s' % (now.strftime('%Y%m%dT%H%M%S'))
        os.mkdir(dirname)

        # make a resources directory and copy the Bootstrap stuff there
        self.copy_resources(dirname)

        # where to find the template files
        templates_dir = os.path.dirname(templates.__file__)
        mylookup = TemplateLookup(directories=[templates_dir])

        # render the simulation steps first
        for stage in self.result.items():
            if 'substages' in stage[1].keys():
                for substage in stage[1]['substages'].items():
                    print 'rendering', stage[0], substage[0]
                    context = {'dirname':dirname, 'renderstage':stage[0],
                      'rendersubstage':substage[0], 'data':self.result}

                    link = os.path.join(context['dirname'], '%s-%s.html' % (stage[0], substage[0]))
                    with open(link, 'w') as f:
                        template = mylookup.get_template('%s-%s.html' % (stage[0], substage[0]))
                        f.write(template.render(**context))

            else:
                # no substages
                print 'rendering', stage[0]
                context = {'mylookup':mylookup, 'dirname':dirname,
                  'renderstage':stage[0],
                  'data':self.result}

                link = os.path.join(context['dirname'], '%s.html' % (stage[0]))
                with open(link, 'w') as f:
                    template = mylookup.get_template('%s.html' % (stage[0]))
                    f.write(template.render(**context))

