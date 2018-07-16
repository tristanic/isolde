/* minimol_io.cpp: atomic model types */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#include "minimol_io.h"

#include <fstream>
#include <sstream>
extern "C" {
#include <string.h>
}


namespace clipper {


// MMDBfile


// local types for referencing between MMDB and minimol
struct SAtom { mmdb::PCAtom db; const MAtom* mm; };
struct SMono { mmdb::PCResidue db; const MMonomer* mm; std::vector<SAtom> data; };
struct SPoly { mmdb::PCChain db; const MPolymer* mm; std::vector<SMono> data; };
struct SModl { mmdb::PCModel db; const MModel* mm; std::vector<SPoly> data; };


/*! The file may be either a PDB or mmCIF file.
  If the spacegroup or cell are not set, they will be taken from the
  file, otherwise the existing values override the ones in the file.
  \param file The filename (or pathname) of the file to read. */
void MMDBfile::read_file( const String& file )
{
  int err = ReadCoorFile( (char *)file.c_str() );
  if (err) Message::message( Message_fatal( "MMDBfile: read_file error: "+file+" : "+String(err) ) );
}

/*! The output file type will be the same as the file read, otherwise PDB.
  \param file The filename (or pathname) of the file to write.
  \param type 0=PDB, 1=CIF, 2=binary, default=same as input file, or PDB. */
void MMDBfile::write_file( const String& file, TYPE type )
{
  const TYPE types[3] = { PDB, CIF, Binary };
  int rtype = GetFileType();
  if ( type == Default && rtype >= 0 && rtype <= 2 ) type = types[ rtype ];
  int err;
  switch ( type ) {
  case Binary:
    err = WriteMMDBF( (char *)file.c_str() ); break;
  case CIF:
    err = WriteCIFASCII( (char *)file.c_str() ); break;
  case PDB:
  default:
    err = WritePDBASCII( (char *)file.c_str() ); break;
  }
  if (err) Message::message( Message_fatal( "MMDBfile: write_file error: "+file+" : "+String(err) ) ); 
}

/*! Import data from the MMDB hierarchy into a MiniMol object. All the
  atoms, residues and chains are imported unless a selection handle is
  supplied, in which case only the selected atoms and their containers
  are imported.

  Each data is tagged with an optional Property<String> (see
  clipper::PropertyManager), called "CID", which describes the
  position in the MMDB hierarchy from which it was taken. This may
  (optionally) be used later to put a modified element back into the
  hierarchy.

  Any atoms which have alternate conformation codes will also be given
  an "AltConf" Property<String>.

  For more details, see clipper::PropertyManager::exists_property(),
  clipper::PropertyManager::get_property().

  \param minimol The minimol to import.
  \param hnd (optional) MMDB selection handle. */
void MMDBfile::import_minimol( MiniMol& minimol, const int hnd )
{
  // clear the model
  minimol = MiniMol( spacegroup(), cell() );

  // make atom, residue, chain selections
  int h_atm = NewSelection();
  int h_res = NewSelection();
  int h_chn = NewSelection();
  int h_mod = NewSelection();
  if ( hnd < 0 ) {
    SelectAtoms( h_atm, 0, 0, ::mmdb::SKEY_NEW );
  } else {
    Select( h_atm, ::mmdb::STYPE_ATOM, hnd, ::mmdb::SKEY_NEW );
  }
  Select( h_res, ::mmdb::STYPE_RESIDUE, h_atm, ::mmdb::SKEY_NEW );
  Select( h_chn, ::mmdb::STYPE_CHAIN,   h_atm, ::mmdb::SKEY_NEW );
  Select( h_mod, ::mmdb::STYPE_MODEL,   h_atm, ::mmdb::SKEY_NEW );

  // now import objects
  char txt[256];
  MModel& mol = minimol.model();
  mmdb::PCModel p_mod = GetModel(1);
  if ( p_mod != NULL ) {
    for ( int c = 0; c < p_mod->GetNumberOfChains(); c++ ) {
      mmdb::PCChain p_chn = p_mod->GetChain(c);
      if ( p_chn != NULL ) if ( p_chn->isInSelection( h_chn ) ) {
	// import the chain
	MPolymer pol;
	for ( int r = 0; r < p_chn->GetNumberOfResidues(); r++ ) {
	  mmdb::PCResidue p_res = p_chn->GetResidue(r);
	  if ( p_chn != NULL ) if ( p_res->isInSelection( h_res ) ) {
	    // import the residue
	    MMonomer mon;
	    for ( int a = 0; a < p_res->GetNumberOfAtoms(); a++ ) {
	      mmdb::PCAtom p_atm = p_res->GetAtom(a);
	      if ( p_atm != NULL ) if ( p_atm->isInSelection( h_atm ) )
	       if ( !p_atm->Ter ) {
		// import the atom
		MAtom atm( Atom::null() );
		atm.set_name( p_atm->GetAtomName(), p_atm->altLoc );
		atm.set_element( p_atm->element );
		if ( p_atm->WhatIsSet & ::mmdb::ASET_Coordinates )
		  atm.set_coord_orth(
		    Coord_orth( p_atm->x, p_atm->y, p_atm->z ) );
		if ( p_atm->WhatIsSet & ::mmdb::ASET_Occupancy )
		  atm.set_occupancy( p_atm->occupancy );
		if ( p_atm->WhatIsSet & ::mmdb::ASET_tempFactor )
		  atm.set_u_iso( Util::b2u( p_atm->tempFactor ) );
		if ( p_atm->WhatIsSet & ::mmdb::ASET_Anis_tFac )
		  atm.set_u_aniso_orth(
		    U_aniso_orth( p_atm->u11, p_atm->u22, p_atm->u33,
				  p_atm->u12, p_atm->u13, p_atm->u23 ) );
		p_atm->GetAtomID( txt );
		atm.set_property("CID",Property<String>(String(txt)));
		if ( p_atm->altLoc[0] != '\0' )
		  atm.set_property("AltConf",
				   Property<String>(String(p_atm->altLoc)));
		mon.insert( atm );  // store the atom
	      }
	    }
	    mon.set_seqnum( p_res->GetSeqNum(), String(p_res->GetInsCode()) );
	    mon.set_type( p_res->GetResName() );
	    p_res->GetResidueID( txt );
	    mon.set_property("CID",Property<String>(String(txt)));
	    pol.insert( mon );  // store the residue
	  }
	}
	pol.set_id( p_chn->GetChainID() );
	p_chn->GetChainID( txt );
	pol.set_property("CID",Property<String>(String(txt)));
	mol.insert( pol );  // store the chain
      }
    }
    p_mod->GetModelID( txt );
    mol.set_property("CID",Property<String>(String(txt)));
  }

  // clean up
  DeleteSelection( h_atm );
  DeleteSelection( h_res );
  DeleteSelection( h_chn );
  DeleteSelection( h_mod );
}


/*! Export data to the MMDB hierarchy from a MiniMol object. All the
  atoms, residues and chains are exported.

  If any MiniMol object has a "CID" Property<String> (see
  clipper::PropertyManager), then the information from that object
  will be used to update the corresponding object in the MMDB
  hierarchy, if it exists. If there is no such entry in the MMDB
  hierarchy, or if no "CID" Property<String> exists, then a new object
  will be created in the MMDB hierarchy. */
void MMDBfile::export_minimol( const MiniMol& minimol )
{

  // export spacegroup/cell
  if ( !minimol.spacegroup().is_null() ) set_spacegroup( minimol.spacegroup() );
  if ( !minimol.cell().is_null() )       set_cell( minimol.cell() );

  // create structure for relationships between Minimol and MMDB
  SModl smod;
  clipper::String cid;

  // fill structure
  smod.mm = &(minimol.model());
  smod.db = NULL;
  smod.data.resize( smod.mm->size() );
  for ( int p = 0; p < smod.data.size(); p++ ) {  // loop over chains
    SPoly& spol = smod.data[p];
    spol.mm = &((*smod.mm)[p]);
    spol.db = NULL;
    spol.data.resize( spol.mm->size() );
    for ( int r = 0; r < spol.data.size(); r++ ) {  // loop over residues
      SMono& smon = spol.data[r];
      smon.mm = &((*spol.mm)[r]);
      smon.db = NULL;
      smon.data.resize( smon.mm->size() );
      for ( int a = 0; a < smon.data.size(); a++ ) {  // loop over atoms
	SAtom& satm = smon.data[a];
	satm.mm = &((*smon.mm)[a]);
	satm.db = NULL;
      }
    }
  }

  // make the MMDB references by CID if present
  if ( smod.mm->exists_property("CID") ) {
    cid = dynamic_cast<const Property<String>&>( smod.mm->get_property("CID") ).value() + "/A/0/A";
    smod.db = GetModel( (char*)cid.c_str() );
  }
  if ( smod.db != NULL )
   for ( int p = 0; p < smod.data.size(); p++ ) {  // loop over chains
    SPoly& spol = smod.data[p];
    if ( spol.mm->exists_property("CID") ) {
      cid = dynamic_cast<const Property<String>&>( spol.mm->get_property("CID") ).value() + "/0/A";
      spol.db = GetChain( (char*)cid.c_str() );
    }
    if ( spol.db != NULL )
     for ( int r = 0; r < spol.data.size(); r++ ) {  // loop over residues
      SMono& smon = spol.data[r];
      if ( smon.mm->exists_property("CID") ) {
	cid = dynamic_cast<const Property<String>&>( smon.mm->get_property("CID") ).value() + "/A";
	smon.db = GetResidue( (char*)cid.c_str() );
      }
      if ( smon.db != NULL )
       for ( int a = 0; a < smon.data.size(); a++ ) {  // loop over atoms
	SAtom& satm = smon.data[a];
	if ( satm.mm->exists_property("CID") ) {
	  cid = dynamic_cast<const Property<String>&>( satm.mm->get_property("CID") ).value();
          satm.db = GetAtom( (char*)cid.c_str() );
	}
       }
     }
   }

  // Now create MMDB objects for anything which is missing
  if ( smod.db == NULL ) {
    smod.db = new mmdb::CModel();
    AddModel( smod.db );
  }
  for ( int p = 0; p < smod.data.size(); p++ ) {  // loop over chains
    SPoly& spol = smod.data[p];
    if ( spol.db == NULL ) {
      spol.db = new mmdb::CChain();
      smod.db->AddChain( spol.db );
    }
    for ( int r = 0; r < spol.data.size(); r++ ) {  // loop over residues
      SMono& smon = spol.data[r];
      if ( smon.db == NULL ) {
	smon.db = new mmdb::CResidue();
	spol.db->AddResidue( smon.db );
      }
      for ( int a = 0; a < smon.data.size(); a++ ) {  // loop over atoms
	SAtom& satm = smon.data[a];
	if ( satm.db == NULL ) {
	  satm.db = new mmdb::CAtom();
	  smon.db->AddAtom( satm.db );
	}
      }
    }
  }

  // now fill in information in mmdb from MiniMol
  for ( int p = 0; p < smod.data.size(); p++ ) {  // loop over chains
    SPoly& spol = smod.data[p];
    spol.db->SetChainID( (char*)spol.mm->id().substr(0,9).c_str() );  // set id
    for ( int r = 0; r < spol.data.size(); r++ ) {  // loop over residues
      SMono& smon = spol.data[r];
      smon.db->seqNum = smon.mm->seqnum();  // set residue info
      smon.db->SetResName( (char*)smon.mm->type().substr(0,19).c_str() );
      int pos = smon.mm->id().find( ":" );
      if ( pos != String::npos )
	strcpy( smon.db->insCode, smon.mm->id().substr(pos+1,9).c_str() );
      for ( int a = 0; a < smon.data.size(); a++ ) {  // loop over atoms
	SAtom& satm = smon.data[a];
        if ( !satm.mm->coord_orth().is_null() ) {       // set atom coord
          satm.db->x = satm.mm->coord_orth().x();
          satm.db->y = satm.mm->coord_orth().y();
          satm.db->z = satm.mm->coord_orth().z();
          satm.db->WhatIsSet |= ::mmdb::ASET_Coordinates;
        }
        if ( !Util::is_nan( satm.mm->occupancy() ) ) {  // set atom occ
          satm.db->occupancy = satm.mm->occupancy();
          satm.db->WhatIsSet |= ::mmdb::ASET_Occupancy;
        }
        if ( !Util::is_nan( satm.mm->u_iso() ) ) {      // set atom u_iso
          satm.db->tempFactor = Util::u2b( satm.mm->u_iso() );
          satm.db->WhatIsSet |= ::mmdb::ASET_tempFactor;
        }
        if ( !satm.mm->u_aniso_orth().is_null() ) {     // set atom u_aniso
          satm.db->u11 = satm.mm->u_aniso_orth()(0,0);
          satm.db->u22 = satm.mm->u_aniso_orth()(1,1);
          satm.db->u33 = satm.mm->u_aniso_orth()(2,2);
          satm.db->u12 = satm.mm->u_aniso_orth()(0,1);
          satm.db->u13 = satm.mm->u_aniso_orth()(0,2);
          satm.db->u23 = satm.mm->u_aniso_orth()(1,2);
          satm.db->WhatIsSet |= ::mmdb::ASET_Anis_tFac;
        }
        if ( satm.mm->id() != "" )       // atom id
          satm.db->SetAtomName( (char*)satm.mm->id().substr(0,19).c_str() );
        if ( satm.mm->element() != "" )  // atom element
          satm.db->SetElementName( (char*)satm.mm->element().substr(0,19).c_str() );
	if ( satm.mm->exists_property("AltConf") ) {    // alt conf code
	  String a = dynamic_cast<const Property<String>&>( satm.mm->get_property("AltConf") ).value();
	  if ( a != "" ) strcpy( satm.db->altLoc, a.substr(0,19).c_str() );
	}
      }
    }
  }
  FinishStructEdit();
}


void SEQfile::read_file( const String& file )
{
  std::ifstream seqfile( file.c_str() );
  std::ostringstream s;
  s << seqfile.rdbuf();
  contents = s.str();
  //non-portable to old Sun
  //contents = std::string(std::istreambuf_iterator<char>(seqfile),
  //           std::istreambuf_iterator<char>());
}


void SEQfile::import_polymer_sequence( MPolymerSequence& target )
{
  MMoleculeSequence mms;
  import_molecule_sequence( mms );
  target = mms[0];
}

void SEQfile::import_molecule_sequence( MMoleculeSequence& target )
{
  std::vector<clipper::String> lines = contents.split("\n");
  clipper::String id, seq = "";
  for ( int l = 0; l < lines.size(); l++ ) {
    clipper::String line = lines[l].trim();
    if ( line[0] == '>' ) {
      if ( seq != "" ) {
	MPolymerSequence s;
	s.set_id( id );
	s.set_sequence( seq );
	target.insert( s );
      }
      id = line.substr(1);
      id = id.trim();
      seq = "";
    } else if ( isalpha(line[0]) ) {
      for ( int i = 0; i < line.length(); i++ )
	if ( isalpha(line[i]) ) seq += toupper(line[i]);
    }
  }
  if ( seq != "" ) {
    MPolymerSequence s;
    s.set_id( id );
    s.set_sequence( seq );
    target.insert( s );
  }
}

} // namespace clipper
