from __future__ import absolute_import

from collections import defaultdict
import os
import pprint
import xlrd

from mako.template import Template

class Base:
    """
    """

    @staticmethod
    def render(context):
        print 'base.render'
        link = os.path.join(context['dirname'], 't2-1-%s.html' % context['result'][0])
        with open(link, 'w') as f:
            t2_1 = Template(filename='t2-1empty.template')
            f.write(t2_1.render(**context))

        return os.path.basename(link)
            
