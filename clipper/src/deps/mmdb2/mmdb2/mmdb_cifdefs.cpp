//  $Id: mmdb_cifdefs.cpp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    21.11.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDBF_Defs <implementation>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Namespace:  mmdb::
//
//      CIF Definitions
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include "mmdb_cifdefs.h"

namespace mmdb  {

  // ------------------------------------------------------------------

  cpstr CIFName ( int NameID, CIF_MODE Mode )  {
  //  Gives CIF name according to CIF Mode.

    switch (Mode)  {

      case CIF_NDB :

        switch (NameID)  {
          case CAT_POLY_SEQ_SCHEME        :
                    return CIFCAT_NDB_POLY_SEQ_SCHEME;
          case TAG_ID_CODE                :
                    return CIFTAG_NDB_PDB_ID_CODE;
          case TAG_CHAIN_ID               :
                    return CIFTAG_NDB_CHAIN_ID;
          case TAG_SEQ_ALIGN_BEG          :
                    return CIFTAG_SEQ_ALIGN_BEG;
          case TAG_SEQ_ALIGN_BEG_INS_CODE :
                    return CIFTAG_NDB_SEQ_ALIGN_BEG_INS_CODE;
          case TAG_SEQ_ALIGN_END          :
                    return CIFTAG_SEQ_ALIGN_END;
          case TAG_SEQ_ALIGN_END_INS_CODE :
                    return CIFTAG_NDB_SEQ_ALIGN_END_INS_CODE;
          case TAG_DB_ACCESSION           :
                    return CIFTAG_NDB_DB_ACCESSION;
          case TAG_DB_ALIGN_BEG           :
                    return CIFTAG_DB_ALIGN_BEG;
          case TAG_DB_ALIGN_BEG_INS_CODE  :
                    return CIFTAG_NDB_DB_ALIGN_BEG_INS_CODE;
          case TAG_DB_ALIGN_END           :
                    return CIFTAG_DB_ALIGN_END;
          case TAG_DB_ALIGN_END_INS_CODE  :
                    return CIFTAG_NDB_DB_ALIGN_END_INS_CODE;
          case TAG_SEQ_CHAIN_ID           :
                    return CIFTAG_ID;
          default : return pstr("ERROR_IN_CIF_NAME_1");
        }

      case CIF_PDBX :

        switch (NameID)  {
          case CAT_POLY_SEQ_SCHEME        :
                    return CIFCAT_PDBX_POLY_SEQ_SCHEME;
          case TAG_ID_CODE                :
                    return CIFTAG_PDBX_PDB_ID_CODE;
          case TAG_CHAIN_ID               :
                    return CIFTAG_PDBX_STRAND_ID;
          case TAG_SEQ_ALIGN_BEG          :
                    return CIFTAG_SEQ_ALIGN_BEG;
          case TAG_SEQ_ALIGN_BEG_INS_CODE :
                    return CIFTAG_PDBX_SEQ_ALIGN_BEG_INS_CODE;
          case TAG_SEQ_ALIGN_END          :
                    return CIFTAG_SEQ_ALIGN_END;
          case TAG_SEQ_ALIGN_END_INS_CODE :
                    return CIFTAG_PDBX_SEQ_ALIGN_END_INS_CODE;
          case TAG_DB_ACCESSION           :
                    return CIFTAG_PDBX_DB_ACCESSION;
          case TAG_DB_ALIGN_BEG           :
                    return CIFTAG_DB_ALIGN_BEG;
          case TAG_DB_ALIGN_BEG_INS_CODE  :
                    return CIFTAG_PDBX_DB_ALIGN_BEG_INS_CODE;
          case TAG_DB_ALIGN_END           :
                    return CIFTAG_DB_ALIGN_END;
          case TAG_DB_ALIGN_END_INS_CODE  :
                    return CIFTAG_PDBX_DB_ALIGN_END_INS_CODE;
          case TAG_SEQ_CHAIN_ID           :
                    return CIFTAG_ASYM_ID;
          default : return pstr("ERROR_IN_CIF_NAME_2");
        }

      default : return pstr("ERROR_IN_CIF_NAME_3");

    }

  }

  cpstr CIFCAT_ATOM_SITE                  = cpstr("_atom_site");
  cpstr CIFCAT_ATOM_SITE_ANISOTROP        = cpstr("_atom_site_anisotrop");
  cpstr CIFCAT_ATOM_SITES                 = cpstr("_atom_sites");
  cpstr CIFCAT_AUDIT_AUTHOR               = cpstr("_audit_author");
  cpstr CIFCAT_CELL                       = cpstr("_cell");
  cpstr CIFCAT_CHEM_COMP                  = cpstr("_chem_comp");
  cpstr CIFCAT_CITATION                   = cpstr("_citation");
  cpstr CIFCAT_DATABASE                   = cpstr("_database");
  cpstr CIFCAT_DATABASE_PDB_CAVEAT        = cpstr("_database_pdb_caveat");
  cpstr CIFCAT_DATABASE_PDB_MATRIX        = cpstr("_database_pdb_matrix");
  cpstr CIFCAT_DATABASE_PDB_REV           = cpstr("_database_pdb_rev");
  cpstr CIFCAT_DATABASE_PDB_TVECT         = cpstr("_database_pdb_tvect");
  cpstr CIFCAT_ENTITY                     = cpstr("_entity");
  cpstr CIFCAT_EXPTL                      = cpstr("_exptl");
  cpstr CIFCAT_NDB_DATABASE_REMARK        = cpstr("_ndb_database_remark");
  cpstr CIFCAT_NDB_NONSTANDARD_LIST       = cpstr("_ndb_nonstandard_list");
  cpstr CIFCAT_NDB_POLY_SEQ_SCHEME        = cpstr("_ndb_poly_seq_scheme");
  cpstr CIFCAT_PDBX_POLY_SEQ_SCHEME       = cpstr("_pdbx_poly_seq_scheme");
  cpstr CIFCAT_REFINE                     = cpstr("_refine");
  cpstr CIFCAT_SPRSDE                     = cpstr("_ndb_database_pdb_obs_spr");
  cpstr CIFCAT_STRUCT                     = cpstr("_struct");
  cpstr CIFCAT_STRUCT_ASYM                = cpstr("_struct_asym");
  cpstr CIFCAT_STRUCT_CONF                = cpstr("_struct_conf");
  cpstr CIFCAT_STRUCT_CONN                = cpstr("_struct_conn");
  cpstr CIFCAT_STRUCT_LINKR               = cpstr("_struct_linkr");
  cpstr CIFCAT_STRUCT_KEYWORDS            = cpstr("_struct_keywords");
  cpstr CIFCAT_STRUCT_NCS_OPER            = cpstr("_struct_ncs_oper");
  cpstr CIFCAT_STRUCT_REF                 = cpstr("_struct_ref");
  cpstr CIFCAT_STRUCT_REF_SEQ             = cpstr("_struct_ref_seq");
  cpstr CIFCAT_STRUCT_REF_SEQ_DIF         = cpstr("_struct_ref_seq_dif");
  cpstr CIFCAT_STRUCT_SHEET               = cpstr("_struct_sheet");
  cpstr CIFCAT_STRUCT_SHEET_RANGE         = cpstr("_struct_sheet_range");
  cpstr CIFCAT_STRUCT_SHEET_ORDER         = cpstr("_struct_sheet_order");
  cpstr CIFCAT_STRUCT_SHEET_HBOND         = cpstr("_struct_sheet_hbond");
  cpstr CIFCAT_SYMMETRY                   = cpstr("_symmetry");
  cpstr CIFCAT_OBSLTE                     = cpstr("_ndb_database_pdb_obs_spr");


  cpstr CIFTAG_ANGLE_ALPHA                   = cpstr("angle_alpha");
  cpstr CIFTAG_ANGLE_BETA                    = cpstr("angle_beta");
  cpstr CIFTAG_ANGLE_GAMMA                   = cpstr("angle_gamma");
  cpstr CIFTAG_ASYM_ID                       = cpstr("asym_id");
  cpstr CIFTAG_ATOM_TYPE_SYMBOL              = cpstr("atom_type_symbol");
  cpstr CIFTAG_AUTH_ASYM_ID                  = cpstr("auth_asym_id");
  cpstr CIFTAG_AUTH_ATOM_ID                  = cpstr("auth_atom_id");
  cpstr CIFTAG_AUTH_COMP_ID                  = cpstr("auth_comp_id");
  cpstr CIFTAG_AUTH_SEQ_ID                   = cpstr("auth_seq_id");
  cpstr CIFTAG_B_ISO_OR_EQUIV                = cpstr("B_iso_or_equiv");
  cpstr CIFTAG_B_ISO_OR_EQUIV_ESD            = cpstr("B_iso_or_equiv_esd");
  cpstr CIFTAG_BEG_LABEL_ASYM_ID             = cpstr("beg_label_asym_id");
  cpstr CIFTAG_BEG_LABEL_COMP_ID             = cpstr("beg_label_comp_id");
  cpstr CIFTAG_BEG_LABEL_SEQ_ID              = cpstr("beg_label_seq_id");
  cpstr CIFTAG_CARTN_X                       = cpstr("cartn_x");
  cpstr CIFTAG_CARTN_X_ESD                   = cpstr("cartn_x_esd");
  cpstr CIFTAG_CARTN_Y                       = cpstr("cartn_y");
  cpstr CIFTAG_CARTN_Y_ESD                   = cpstr("cartn_y_esd");
  cpstr CIFTAG_CARTN_Z                       = cpstr("cartn_z");
  cpstr CIFTAG_CARTN_Z_ESD                   = cpstr("cartn_z_esd");
  cpstr CIFTAG_PDBX_FORMAL_CHARGE            = cpstr("pdbx_formal_charge");
  cpstr CIFTAG_CODE                          = cpstr("code");
  cpstr CIFTAG_CODE_NDB                      = cpstr("code_NDB");
  cpstr CIFTAG_CODE_PDB                      = cpstr("code_PDB");
  cpstr CIFTAG_CONF_TYPE_ID                  = cpstr("conf_type_id");
  cpstr CIFTAG_CONN_TYPE_ID                  = cpstr("conn_type_id");
  cpstr CIFTAG_DATE                          = cpstr("date");
  cpstr CIFTAG_DATE_ORIGINAL                 = cpstr("date_original");
  cpstr CIFTAG_DB_ALIGN_BEG                  = cpstr("db_align_beg");
  cpstr CIFTAG_DB_ALIGN_END                  = cpstr("db_align_end");
  cpstr CIFTAG_DB_CODE                       = cpstr("db_code");
  cpstr CIFTAG_DB_MON_ID                     = cpstr("db_mon_id");
  cpstr CIFTAG_DB_NAME                       = cpstr("db_name");
  cpstr CIFTAG_DETAILS                       = cpstr("details");
  cpstr CIFTAG_END_LABEL_ASYM_ID             = cpstr("end_label_asym_id");
  cpstr CIFTAG_END_LABEL_COMP_ID             = cpstr("end_label_comp_id");
  cpstr CIFTAG_END_LABEL_SEQ_ID              = cpstr("end_label_seq_id");
  cpstr CIFTAG_ENTITY_ID                     = cpstr("entity_id");
  cpstr CIFTAG_ENTRY_ID                      = cpstr("entry_id");
  cpstr CIFTAG_FORMULA                       = cpstr("formula");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX11         = cpstr("fract_transf_matrix[1][1]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX12         = cpstr("fract_transf_matrix[1][2]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX13         = cpstr("fract_transf_matrix[1][3]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX21         = cpstr("fract_transf_matrix[2][1]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX22         = cpstr("fract_transf_matrix[2][2]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX23         = cpstr("fract_transf_matrix[2][3]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX31         = cpstr("fract_transf_matrix[3][1]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX32         = cpstr("fract_transf_matrix[3][2]");
  cpstr CIFTAG_FRACT_TRANSF_MATRIX33         = cpstr("fract_transf_matrix[3][3]");
  cpstr CIFTAG_FRACT_TRANSF_VECTOR1          = cpstr("fract_transf_vector[1]");
  cpstr CIFTAG_FRACT_TRANSF_VECTOR2          = cpstr("fract_transf_vector[2]");
  cpstr CIFTAG_FRACT_TRANSF_VECTOR3          = cpstr("fract_transf_vector[3]");
  cpstr CIFTAG_GROUP_PDB                     = cpstr("group_PDB" );
  cpstr CIFTAG_ID                            = cpstr("id");
  cpstr CIFTAG_INS_CODE                      = cpstr("ins_code");
  cpstr CIFTAG_LABEL_ALT_ID                  = cpstr("label_alt_id");
  cpstr CIFTAG_LABEL_ATOM_ID                 = cpstr("label_atom_id");
  cpstr CIFTAG_LABEL_ASYM_ID                 = cpstr("label_asym_id");
  cpstr CIFTAG_LABEL_COMP_ID                 = cpstr("label_comp_id");
  cpstr CIFTAG_LABEL_ENTITY_ID               = cpstr("label_entity_id");
  cpstr CIFTAG_LABEL_SEQ_ID                  = cpstr("label_seq_id");
  cpstr CIFTAG_LENGTH_A                      = cpstr("length_a");
  cpstr CIFTAG_LENGTH_B                      = cpstr("length_b");
  cpstr CIFTAG_LENGTH_C                      = cpstr("length_c");
  cpstr CIFTAG_LS_D_RES_HIGH                 = cpstr("ls_d_res_high");
  cpstr CIFTAG_MATRIX11                      = cpstr("matrix[1][1]");
  cpstr CIFTAG_MATRIX12                      = cpstr("matrix[1][2]");
  cpstr CIFTAG_MATRIX13                      = cpstr("matrix[1][3]");
  cpstr CIFTAG_MATRIX21                      = cpstr("matrix[2][1]");
  cpstr CIFTAG_MATRIX22                      = cpstr("matrix[2][2]");
  cpstr CIFTAG_MATRIX23                      = cpstr("matrix[2][3]");
  cpstr CIFTAG_MATRIX31                      = cpstr("matrix[3][1]");
  cpstr CIFTAG_MATRIX32                      = cpstr("matrix[3][2]");
  cpstr CIFTAG_MATRIX33                      = cpstr("matrix[3][3]");
  cpstr CIFTAG_METHOD                        = cpstr("method");
  cpstr CIFTAG_MOD_TYPE                      = cpstr("mod_type");
  cpstr CIFTAG_MON_ID                        = cpstr("mon_id");
  cpstr CIFTAG_NAME                          = cpstr("name");
  cpstr CIFTAG_NDB_BEG_LABEL_INS_CODE_PDB    = cpstr("ndb_beg_label_ins_code_pdb");
  cpstr CIFTAG_NDB_CHAIN_ID                  = cpstr("ndb_chain_id");
  cpstr CIFTAG_NDB_COMPONENT_NO              = cpstr("ndb_component_no");
  cpstr CIFTAG_NDB_DESCRIPTOR                = cpstr("ndb_descriptor");
  cpstr CIFTAG_NDB_DB_ACCESSION              = cpstr("ndb_db_accession");
  cpstr CIFTAG_NDB_DB_ALIGN_BEG_INS_CODE     = cpstr("ndb_db_align_beg_ins_code");
  cpstr CIFTAG_NDB_DB_ALIGN_END_INS_CODE     = cpstr("ndb_db_align_end_ins_code");
  cpstr CIFTAG_NDB_END_LABEL_INS_CODE_PDB    = cpstr("ndb_end_label_ins_code_pdb");
  //cpstr CIFTAG_NDB_INS_CODE                =   cpstr("ndb_ins_code");
  cpstr CIFTAG_PDBX_PDB_INS_CODE             = cpstr("pdbx_PDB_ins_code");
  cpstr CIFTAG_NDB_HELIX_CLASS_PDB           = cpstr("ndb_helix_class_pdb");
  cpstr CIFTAG_NDB_KEYWORDS                  = cpstr("ndb_keywords");
  cpstr CIFTAG_NDB_LABEL_ALT_ID              = cpstr("ndb_label_alt_id");
  cpstr CIFTAG_NDB_LABEL_ATOM_ID             = cpstr("ndb_label_atom_id");
  cpstr CIFTAG_NDB_LABEL_ASYM_ID             = cpstr("ndb_label_asym_id");
  cpstr CIFTAG_NDB_LABEL_COMP_ID             = cpstr("ndb_label_comp_id");
  cpstr CIFTAG_NDB_LABEL_INS_CODE            = cpstr("ndb_label_ins_code");
  cpstr CIFTAG_NDB_LABEL_SEQ_NUM             = cpstr("ndb_label_seq_num");
  cpstr CIFTAG_NDB_LENGTH                    = cpstr("ndb_length");
  cpstr CIFTAG_NDB_MODEL                     = cpstr("ndb_model");
  cpstr CIFTAG_NDB_PDB_CHAIN_ID              = cpstr("ndb_pdb_chain_id");
  cpstr CIFTAG_NDB_PDB_ID                    = cpstr("ndb_pdb_id");
  cpstr CIFTAG_NDB_PDB_ID_CODE               = cpstr("ndb_pdb_id_code");
  cpstr CIFTAG_NDB_PDB_INS_CODE              = cpstr("ndb_pdb_ins_code");
  cpstr CIFTAG_NDB_PTNR1_LABEL_INS_CODE      = cpstr("ndb_ptnr1_label_ins_code");
  cpstr CIFTAG_NDB_PTNR1_STANDARD_COMP_ID    = cpstr("ndb_ptnr1_standard_comp_id");
  cpstr CIFTAG_NDB_RANGE_1_BEG_LABEL_COMP_ID = cpstr("ndb_range_1_beg_label_comp_id");
  cpstr CIFTAG_NDB_RANGE_1_BEG_LABEL_ASYM_ID = cpstr("ndb_range_1_beg_label_asym_id");
  cpstr CIFTAG_NDB_RANGE_1_BEG_LABEL_INS_CODE= cpstr("ndb_range_1_beg_label_ins_code");
  cpstr CIFTAG_NDB_RANGE_1_END_LABEL_COMP_ID = cpstr("ndb_range_1_end_label_comp_id");
  cpstr CIFTAG_NDB_RANGE_1_END_LABEL_ASYM_ID = cpstr("ndb_range_1_end_label_asym_id");
  cpstr CIFTAG_NDB_RANGE_1_END_LABEL_INS_CODE= cpstr("ndb_range_1_end_label_ins_code");
  cpstr CIFTAG_NDB_SEQ_ALIGN_BEG             = cpstr("ndb_seq_align_beg");
  cpstr CIFTAG_NDB_SEQ_ALIGN_BEG_INS_CODE    = cpstr("ndb_seq_align_beg_ins_code");
  cpstr CIFTAG_NDB_SEQ_ALIGN_END             = cpstr("ndb_seq_align_end");
  cpstr CIFTAG_NDB_SEQ_ALIGN_END_INS_CODE    = cpstr("ndb_seq_align_end_ins_code");
  cpstr CIFTAG_NDB_SEQ_DB_NAME               = cpstr("ndb_seq_db_name");
  cpstr CIFTAG_NDB_SEQ_DB_ACCESSION_CODE     = cpstr("ndb_seq_db_accession_code");
  cpstr CIFTAG_NDB_SEQ_DB_SEQ_NUM            = cpstr("ndb_seq_db_seq_num");
  cpstr CIFTAG_NDB_SYNONYMS                  = cpstr("ndb_synonyms");
  cpstr CIFTAG_NUM                           = cpstr("num");
  cpstr CIFTAG_NUMBER_ATOMS_NH               = cpstr("number_atoms_nh");
  cpstr CIFTAG_NUMBER_STRANDS                = cpstr("number_strands");
  cpstr CIFTAG_OCCUPANCY                     = cpstr("occupancy");
  cpstr CIFTAG_OCCUPANCY_ESD                 = cpstr("occupancy_esd");
  cpstr CIFTAG_ORIGX11                       = cpstr("origx[1][1]");
  cpstr CIFTAG_ORIGX12                       = cpstr("origx[1][2]");
  cpstr CIFTAG_ORIGX13                       = cpstr("origx[1][3]");
  cpstr CIFTAG_ORIGX21                       = cpstr("origx[2][1]");
  cpstr CIFTAG_ORIGX22                       = cpstr("origx[2][2]");
  cpstr CIFTAG_ORIGX23                       = cpstr("origx[2][3]");
  cpstr CIFTAG_ORIGX31                       = cpstr("origx[3][1]");
  cpstr CIFTAG_ORIGX32                       = cpstr("origx[3][2]");
  cpstr CIFTAG_ORIGX33                       = cpstr("origx[3][3]");
  cpstr CIFTAG_ORIGX_VECTOR1                 = cpstr("origx_vector[1]");
  cpstr CIFTAG_ORIGX_VECTOR2                 = cpstr("origx_vector[2]");
  cpstr CIFTAG_ORIGX_VECTOR3                 = cpstr("origx_vector[3]");
  cpstr CIFTAG_PDB_ID                        = cpstr("pdb_id");
  cpstr CIFTAG_PDB_MON_ID                    = cpstr("pdb_mon_id");
  cpstr CIFTAG_PDB_STRAND_ID                 = cpstr("pdb_strand_id");
  cpstr CIFTAG_PDBX_DB_ACCESSION             = cpstr("pdbx_db_accession");
  cpstr CIFTAG_PDBX_DB_ALIGN_BEG_INS_CODE    = cpstr("pdbx_db_align_beg_ins_code");
  cpstr CIFTAG_PDBX_DB_ALIGN_END_INS_CODE    = cpstr("pdbx_db_align_end_ins_code");
  cpstr CIFTAG_PDBX_PDB_ID_CODE              = cpstr("pdbx_PDB_id_code");
  cpstr CIFTAG_PDBX_PDB_MODEL_NUM            = cpstr("pdbx_PDB_model_num");
  cpstr CIFTAG_PDBX_STRAND_ID                = cpstr("pdbx_strand_id");
  cpstr CIFTAG_RANGE_1_BEG_LABEL_ATOM_ID     = cpstr("range_1_beg_label_atom_id");
  cpstr CIFTAG_RANGE_1_BEG_LABEL_SEQ_ID      = cpstr("range_1_beg_label_seq_id");
  cpstr CIFTAG_RANGE_1_END_LABEL_ATOM_ID     = cpstr("range_1_end_label_atom_id");
  cpstr CIFTAG_RANGE_1_END_LABEL_SEQ_ID      = cpstr("range_1_end_label_seq_id");
  cpstr CIFTAG_RANGE_ID_1                    = cpstr("range_id_1");
  cpstr CIFTAG_RANGE_ID_2                    = cpstr("range_id_2");
  cpstr CIFTAG_RCSB_RECORD_REVISED_1         = cpstr("rcsb_record_revised_1");
  cpstr CIFTAG_RCSB_RECORD_REVISED_2         = cpstr("rcsb_record_revised_2");
  cpstr CIFTAG_RCSB_RECORD_REVISED_3         = cpstr("rcsb_record_revised_3");
  cpstr CIFTAG_RCSB_RECORD_REVISED_4         = cpstr("rcsb_record_revised_4");
  cpstr CIFTAG_PDBX_SEQ_ALIGN_BEG_INS_CODE   = cpstr("pdbx_seq_align_beg_ins_code");
  cpstr CIFTAG_PDBX_SEQ_ALIGN_END_INS_CODE   = cpstr("pdbx_seq_align_end_ins_code");
  cpstr CIFTAG_PTNR1_LABEL_ASYM_ID           = cpstr("ptnr1_label_asym_id");
  cpstr CIFTAG_PTNR1_LABEL_COMP_ID           = cpstr("ptnr1_label_comp_id");
  cpstr CIFTAG_PTNR1_LABEL_SEQ_ID            = cpstr("ptnr1_label_seq_id");
  cpstr CIFTAG_REF_ID                        = cpstr("ref_id");
  cpstr CIFTAG_REPLACES                      = cpstr("replaces");
  cpstr CIFTAG_REPLACE_PDB_ID                = cpstr("replace_pdb_id");
  cpstr CIFTAG_SEGMENT_ID                    = cpstr("segment_id");
  cpstr CIFTAG_SEQ_ALIGN_BEG                 = cpstr("seq_align_beg");
  cpstr CIFTAG_SEQ_ALIGN_END                 = cpstr("seq_align_end");
  cpstr CIFTAG_SEQ_NUM                       = cpstr("seq_num");
  cpstr CIFTAG_SENSE                         = cpstr("sense");
  cpstr CIFTAG_SHEET_ID                      = cpstr("sheet_id");
  cpstr CIFTAG_SOURCE                        = cpstr("source");
  cpstr CIFTAG_SPACE_GROUP_NAME_H_M          = cpstr("space_group_name_H-M");
  cpstr CIFTAG_TEXT                          = cpstr("text");
  cpstr CIFTAG_TITLE                         = cpstr("title");
  cpstr CIFTAG_TYPE                          = cpstr("type");
  cpstr CIFTAG_TYPE_SYMBOL                   = cpstr("type_symbol");
  cpstr CIFTAG_VECTOR1                       = cpstr("vector[1]");
  cpstr CIFTAG_VECTOR2                       = cpstr("vector[2]");
  cpstr CIFTAG_VECTOR3                       = cpstr("vector[3]");
  cpstr CIFTAG_U11                           = cpstr("u[1][1]");
  cpstr CIFTAG_U11_ESD                       = cpstr("u[1][1]_esd");
  cpstr CIFTAG_U12                           = cpstr("u[1][2]");
  cpstr CIFTAG_U12_ESD                       = cpstr("u[1][2]_esd");
  cpstr CIFTAG_U13                           = cpstr("u[1][3]");
  cpstr CIFTAG_U13_ESD                       = cpstr("u[1][3]_esd");
  cpstr CIFTAG_U22                           = cpstr("u[2][2]");
  cpstr CIFTAG_U22_ESD                       = cpstr("u[2][2]_esd");
  cpstr CIFTAG_U23                           = cpstr("u[2][3]");
  cpstr CIFTAG_U23_ESD                       = cpstr("u[2][3]_esd");
  cpstr CIFTAG_U33                           = cpstr("u[3][3]");
  cpstr CIFTAG_U33_ESD                       = cpstr("u[3][3]_esd");
  cpstr CIFTAG_Z_PDB                         = cpstr("z_pdb");

  cpstr CIFTAG_CONN_PTNR1_AUTH_ATOM_ID       = cpstr("ptnr1_auth_atom_id");
  cpstr CIFTAG_CONN_PDBX_PTNR1_AUTH_ALT_ID   = cpstr("pdbx_ptnr1_auth_alt_id");
  cpstr CIFTAG_CONN_PTNR1_AUTH_COMP_ID       = cpstr("ptnr1_auth_comp_id");
  cpstr CIFTAG_CONN_PTNR1_AUTH_ASYM_ID       = cpstr("ptnr1_auth_asym_id");
  cpstr CIFTAG_CONN_PTNR1_AUTH_SEQ_ID        = cpstr("ptnr1_auth_seq_id");
  cpstr CIFTAG_CONN_PDBX_PTNR1_PDB_INS_CODE  = cpstr("pdbx_ptnr1_PDB_ins_code");
  cpstr CIFTAG_CONN_DIST                     = cpstr("link_dist");
  cpstr CIFTAG_CONN_PTNR2_AUTH_ATOM_ID       = cpstr("ptnr2_auth_atom_id");
  cpstr CIFTAG_CONN_PDBX_PTNR2_AUTH_ALT_ID   = cpstr("pdbx_ptnr2_auth_alt_id");
  cpstr CIFTAG_CONN_PTNR2_AUTH_COMP_ID       = cpstr("ptnr2_auth_comp_id");
  cpstr CIFTAG_CONN_PTNR2_AUTH_ASYM_ID       = cpstr("ptnr2_auth_asym_id");
  cpstr CIFTAG_CONN_PTNR2_AUTH_SEQ_ID        = cpstr("ptnr2_auth_seq_id");
  cpstr CIFTAG_CONN_PDBX_PTNR2_PDB_INS_CODE  = cpstr("pdbx_ptnr2_PDB_ins_code");
  cpstr CIFTAG_CONN_PTNR1_SYMMETRY           = cpstr("ptnr1_symmetry");
  cpstr CIFTAG_CONN_PTNR2_SYMMETRY           = cpstr("ptnr2_symmetry");
  cpstr CIFTAG_CONN_NAME                     = cpstr("link_name");

}  // namespace mmdb
