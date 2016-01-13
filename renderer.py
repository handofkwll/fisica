"""This module contains classes and methods used to oversee the rendering
of the web log.
"""

from __future__ import absolute_import

import collections
import os
import shutil

from extern import mako
from mako.template import Template
from mako.lookup import TemplateLookup

import templates
import templates.resources as resources


class Renderer(object):
    """Class to oversee rendering of results into html.
    """

    def __init__(self, result):
        """Constructor.

        Keyword arguments:
          result - structure holding the results of the simulation run.
        """
        self.result = result

    def copy_resources(self, dirname):
        """Method to copy web resources into the target directory.

        Keyword arguments:
          dirname -- target directory.
        """
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

    def run(self, prefix=''):
        """Method called to do the work.

        Keyword arguments:
        prefix -- string to insert into name of directory holding the
                  weblog. Used to differentiate between simulator and
                  reduction results when using the same code to write 
                  the weblogs. Default ''. 
        """

        # make a directory to contain this result
        runtime = self.result['runtime']
        if prefix=='':
            dirname = 'fisica-%s' % runtime
        else:
            dirname = 'fisica-%s-%s' % (prefix, runtime)
        os.mkdir(dirname)
        print 'Weblog:'
        print '  weblog written to %s' % dirname

        # make a resources directory and copy the Bootstrap stuff there
        self.copy_resources(dirname)

        # where to find the template files
        templates_dir = os.path.dirname(templates.__file__)
        mylookup = TemplateLookup(directories=[templates_dir])

        # render each entry in the result structure, except for some 
        # to be ignored
        ignore = ['runtime', 'timeline'] 
        for stage in self.result.items():

            if stage[0] in ignore:
                pass

            elif isinstance(stage[1], dict) and 'substages' in stage[1].keys():
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

