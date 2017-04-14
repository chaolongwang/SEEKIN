/*
 *  Copyright (C) 2015-2015  Genome institute of Singapore (GIS),
 *                           Jinzhuang Dou and Chaolong Wang
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include "extfunc.h"

using namespace std;
using namespace arma;

/*************************************************
Function:
	Main function 
Usage:
    @see  display_usage()

*************************************************/

extern int modelAF ( int arc, char ** argv );
extern int getAF (int arc, char ** argv);
extern int kinship (int argc, char ** argv);

static void pipeline ( int argc, char ** argv );
bool seekin_display_usage();


int main (int argc, char ** argv ){

	fprintf ( stderr, "\nVersion 1.00: released on 2017-3-10 \n");
	argc--;
	argv++;
	if ( argc == 0 )
	{
		seekin_display_usage ();
		return 0;
	}

	if ( strcmp ( "modelAF", argv[0] ) == 0 )
	{
		modelAF ( argc, argv );
	}
	else if(strcmp ( "getAF", argv[0] ) == 0 ){
		getAF ( argc, argv );

	}
	else if(strcmp ( "kinship", argv[0] ) == 0 ){
		kinship (argc, argv );
	}
	else
	{
		seekin_display_usage ();
	}
	return 1;
}



static void pipeline ( int argc, char ** argv )
{
	fprintf ( stderr, "    no function available now for this option\n");
}

bool seekin_display_usage ()
{
	fprintf ( stderr, "\nUsage: seekin <command> [option]\n" );
	fprintf ( stderr, "    modelAF	Model allele frequencies as a linear function of PCs\n" );
	fprintf ( stderr, "    getAF	Estimate individual-specific allele frequencies\n");
	fprintf ( stderr, "    kinship	Estimate kinship and inbreeding coefficients\n");
}


