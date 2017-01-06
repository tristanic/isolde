#
#  Copyright 2016 Jon Agirre & The University of York
#  Distributed under the terms of the LGPL (www.fsf.org)
#
#  Package containing callback functions for dealing with logs & XML
#  To do: add other callback functions

from datetime import datetime
from lxml import etree

CHAR_COUNTER = 0


def setup_callback () :

    global CHAR_COUNTER
    CHAR_COUNTER = 0


def interactive_flush ( log_string="", xml_tree=None  ) :

    global CHAR_COUNTER
    print log_string[CHAR_COUNTER:len(log_string)]
    CHAR_COUNTER = len(log_string)


def offline_flush ( log_string="", xml_tree=None ) :

    current_date = datetime.now().strftime('%Y/%m/%d %H:%M.%S')
    
    with open("log.txt", "w") as log_file:
        log_file.write("{0}\n{1}".format(current_date, log_string))
    
    with open("program.xml", "w") as xml_file:
        string_tree = etree.tostring( xml_tree, pretty_print=True )
        xml_file.write(string_tree)
