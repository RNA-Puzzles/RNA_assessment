#  msgs.py
#  
#  Copyright 2017 Chichau Miau <zmiao@ebi.ac.uk>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  a util function to show the messages

import sys

def show( prefix, message, new_line=True, back=False ):
    new_line_txt = new_line and "\n" or ""
    back_txt = back and "\b" * 50 or ""
    
    sys.stderr.write( "%s%-10s >> %s%s" %(back_txt, prefix[:10], message, new_line_txt) )
    
    if( back ):
        sys.stderr.flush()
    
    if( prefix == "FATAL" ):
        sys.stderr.write( "Abrupt termination!" )
        quit()
