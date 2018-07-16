//  $Id: mmdb_mmcif_.h $
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
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_MMCIF <interface>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::mmcif::Category ( mmCIF category    )
//       ~~~~~~~~~  mmdb::mmcif::Struct   ( mmCIF structure   )
//                  mmdb::mmcif::Loop     ( mmCIF loop        )
//                  mmdb::mmcif::Data     ( mmCIF data block  )
//                  mmdb::mmcif::File     ( mmCIF file        )
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//


#ifndef __MMDB_MMCIF__
#define __MMDB_MMCIF__


#include "mmdb_io_stream.h"
#include "imex.h"

namespace mmdb  {

  namespace mmcif  {


    //  ======================  Category  ==========================

    enum MMCIF_ITEM  {
      MMCIF_Category = 0,
      MMCIF_Struct   = 1,
      MMCIF_Loop     = 2,
      MMCIF_Data     = 3
    };

    DefineClass(Category);
    DefineStreamFunctions(Category);

    /// \brief mmcif::Category is a base class for mmcif::Struct and
    ///        mmcif::Loop, implementations of mmCIF's "structure" and
    ///        "loop" categories.
    /*!
    This class is not instantiated independently in any applications,
    however, it provides a few public functions which work for
    both mmcif::Struct and mmcif::Loop.

    All data in mmCIF hierarchy is addressed using construct
    "category.tag" plus row number (>=0) for loops. Category names
    should always start from underscore, while tags normally start
    with a letter, e.g. "_barrel.id".

    See general principles of working with mmCIF files and mmCIF
    hierarchies in Section \"\ref mmcif_handler\".
    */

    class MMDB_IMEX Category : public io::Stream  {

      friend class Data;

      public :

        /// \brief Basic constructor.
        Category ();

        /// \brief Constructor that assigns category name.
        /// \param[in] N category name (must start with underscore).
        Category ( cpstr N );

        /// \brief Constructor for MMDB data streaming functions.
        Category ( io::RPStream Object );

        /// \brief Destructor.
        ~Category();

        /// \brief Returns category name.
        /// \return NULL if name was not set
        /// \return pointer to character string if name was set
        inline pstr   GetCategoryName() { return name; }

        /// \brief Sets category name.
        /// \param N new category name
        void   SetCategoryName ( cpstr N );

        /// \brief Returns category type.
        /// This function may be used when retrieving categories
        /// (structures and loops) from data blocks (mmcif::Data).
        /// \return MMCIF_Category for mmcif::Category
        /// \return MMCIF_Struct   for mmcif::Struct
        /// \return MMCIF_Loop     for mmcif::Loop
        virtual MMCIF_ITEM GetCategoryID() { return MMCIF_Category; }

        /// \brief Virtual function for writing category's content
        ///        into mmCIF file.
        /// Default implementation does nothing.
        virtual void WriteMMCIF ( io::RFile ) {}

        /// \brief Virtual function for optimizig data structures.
        /// Optimized data structures take less RAM and their indexes
        /// are sorted for quicker access. Sorting is done automatically
        /// as new data is added to the category. If the
        /// category is edited (fields/data removed), it may need
        /// optimization and re-sorting for efficiency.\n\n
        /// The sorting preserves the order of actual appearance of
        /// tags in mmCIF file. If a category is created
        /// programmatically, the order of tags in mmCIF file will be
        /// the same as order of adding them to the category.
        virtual void Optimize();

        /// \brief Sorts category's data for quicker access.
        /// The sorting is essentially re-indexing of data for quicker
        /// access. It does not change the order of data in both mmCIF
        /// hierarchy and mmCIF file. E.g., if tag "serial_no" was 2nd
        /// one in given category before sorting, it will remain on 2nd
        /// place after it, therefore no change in tag number passed
        /// to functions in mmcif::Struct, mmcif::Loop and mmcif::Data.
        void  Sort();

        /// \brief Returns serial number of a tag in the category.
        /// \param[in]  ttag tag (or name of a field in category)
        /// \return \b >=0 : the tag is in given position
        /// \return \b <0 : the tag was not found, but it could be
        ///         inserted before tag with (-rc-1)th index, where
        ///         'rc' is the return.
        int   GetTagNo ( cpstr ttag );

        /// \brief Adds a tag to the category.
        /// Adding a tag in mmcif::Category does not reserve any
        /// placeholder for the corresponding value. All tags get
        /// automatically sorted (reindexed) for quicker access.
        /// Tags will appear in mmCIF file in order of their addition
        /// to the category.
        /// \param[in]  ttag tag to be added.
        /// \return \b >=0 the tag is already in the category, and return
        ///             is its serial number. No changes to the category
        ///             is done
        /// \return \b <0  the tag was added to the list of tags, and
        ///             return is minus total number of tags in the
        ///             category.
        int   AddTag ( cpstr ttag );

        /// \brief Returns the total number of tags in the category
        int   GetNofTags() { return nTags; }

        /// \brief Returns tag with the specified serial number.
        /// The tags are enumerated as 0..GetNofTags()-1.
        /// \param tagNo tag's serial number
        /// \return \b NULL: tagNo is outside the range
        ///                  of 0..GetNofTags()-1
        /// \return \b not \b NULL: tag in tagNo'th position
        pstr  GetTag ( int tagNo );  // 0..nTags-1

        /// \brief Prints list of tags to stdout.
        /// Both sorted and unsorted tags are printed to standard
        /// output. This function may be used for debugging.
        void  PrintTags();

        /// \brief Returns true if all tags from the list are found
        ///        in the category.
        /// The size of the list of tags may be less than the number
        /// of tags in the category, and order of tags is disregarded.
        /// \param[in] tagList  list of tags to be checked for presence
        ///                 in the category. The list must end with NULL
        ///                 pointer, or your program will crash.
        /// \return \b true  if all tags from the list were found in the
        ///               category
        /// \return \b false if one or more tags from the list were not
        ///               found in the category.
        ///
        /// Example:
        /// \code
        ///  cpstr tagList[] = {"id","type","date",NULL};
        ///  mmcif::Struct cifStruct;
        ///  if (cifStruct.CheckTags(tagList))
        ///    printf ( " all tags are found in category %s\n",
        ///             cifStruct.GetCategoryName() );
        /// \endcode
        /// This function is useful for finding categories in
        /// "dirty cifs", where category name is not given.
        bool CheckTags ( cpstr * tagList );

        /// \brief Deep copy of categories.
        /// Deep copy duplicates all data and memory allocations,
        /// producing a genuine clone of the original. Only deep copy
        /// should be used for copying MMDB objects, a mere assignment
        /// operator will fail you.
        /// \param[in] Category a pointer to mmcif::Category, the content of
        ///                 which is copied into 'this' category.
        virtual void Copy ( PCategory Category );

        /// \brief MMDB stream writer.
        void  write ( io::RFile f );

        /// \brief MMDB stream reader.
        void  read  ( io::RFile f );

      protected:
        int      nTags;
        pstr     name;
        psvector tag;
        ivector  index;
        int      nAllocTags;

        void          InitCategory    ();
        virtual void  FreeMemory      ();
        void          ExpandTags      ( int nTagsNew );
        void          PutCategoryName ( cpstr newName );

    };



    //  ======================  Struct  ============================

    DefineClass(Struct);
    DefineStreamFunctions(Struct);

    /// \brief Constants used to specify mmCIF's \"data not given\" and
    /// \"data not available\" data types.
    extern const int CIF_NODATA_DOT;
    extern const int CIF_NODATA_QUESTION;
    extern cpstr     CIF_NODATA_DOT_FIELD;
    extern cpstr     CIF_NODATA_QUESTION_FIELD;

    /// \brief mmcif::Struct represents mmCIF's \"structure\" category,
    ///        where data follows directly the corresponding tag.
    /*!
    mmCIF's \"structure\" category has the following form:
    \code
    _structure_name.tag0   value0
    _structure_name.tag1   value1
    ...........
    _structure_name.tagN   valueN
    \endcode
    mmcif::Struct represents this construct by keeping category name
    (\"_structure_name\") and associated lists of tags
    (\"tag0,tag1...tagN\") and their values (\"value0,value1...valueN\").

    The structure is created automatically when an mmCIF file is read,
    or it may be created programatically and then pushed into file.

    Access to data is provided via tags. Internally, all values are kept
    as character fields, and it is only on the retrieval stage that they
    are converted to other data types (integers, floats or strings).
    If conversion is not possible, an error code is returned by the
    corresponding functions, which should be checked by the application.

    See general principles of working with mmCIF files and mmCIF
    hierarchies, as well as some code samples, in Section
    \"\ref mmcif_handler\".
    */

    class MMDB_IMEX Struct : public Category  {

      public :

        /// \brief Basic constructor
        Struct ();

        /// \brief Constructor that assigns structure name.
        /// \param[in] N structure name (must start with underscore).
        Struct ( cpstr N );

        /// \brief Constructor for MMDB data streaming functions
        Struct ( io::RPStream Object );

        /// \brief Destructor
        ~Struct();

        /// \brief Adds field to the structure.
        /// \param[in] F field value
        /// \param[in] T tag name
        /// \param[in] Concatenate flag to concatenate existing field
        ///            with the value of \b F. If tag \b T is already in
        ///            the structure and \b Concatenate=true, then
        ///            value of \b F is appended to the existing field.
        ///            Otherwise, the field is replaced with the value
        ///            of \b F
        void AddField ( cpstr F, cpstr T, bool Concatenate=false );

        /// \brief Returns category type \b MMCIF_Struct.
        MMCIF_ITEM  GetCategoryID()  { return MMCIF_Struct; }

        /// \brief Optimizes structure for RAM and data access speed.
        /// Optimized data structures take less RAM and their indexes
        /// are sorted for quicker access. Sorting is done automatically
        /// as new data is added to the category. If the structure
        /// is edited (fields/data removed), it may need
        /// optimization and re-sorting for efficiency.\n\n
        /// The sorting preserves the order of actual appearance of
        /// tags in mmCIF file. If a structure is created
        /// programmatically, the order of tags in mmCIF file will be
        /// the same as order of adding them to the structure.
        void Optimize();

        /// \brief Returns value of field corresponding to tag in the
        ///        specified position.
        /// Tag positions are defined by the order of their appearance in
        /// mmCIF file (if structure was read from a file), or by the
        /// order of their addition to the structure (if structure was
        /// created programmatically). Tags are numbered as
        /// 0...GetNofTags()-1.
        /// \param[in] tagNo tag number (position in the structure)
        /// \return \b NULL: tag does not exist
        /// \return \b CIF_NODATA_DOT_FIELD the field contains
        ///            \"data not given\" value
        /// \return \b CIF_NODATA_QUESTION_FIELD the field contains
        ///            \"data not available\" value
        /// \return \b not \b NULL: string value of the field
        pstr GetField ( int tagNo );  // 0..nTags-1

        /// \brief Fetches value, corresponding to the given tag, as
        ///        a string
        /// \param[out] S pointer to string, which will point to newly
        ///               allocated character string, containing value
        ///               associated with tag \b TName. If tag or value
        ///               is not found, or if value corresponds to
        ///               mmCIF's \"data not given\" or
        ///               \"data not available\", \b S returns NULL.
        /// \param[in] TName character string with tag name
        /// \param[in] Remove flag to remove the tag and its value from
        ///               structure after it is read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_NoField: value is not found
        /// \return \b CIFRC_Ok: success. If \b S returns NULL, then
        ///                the value corresponds to either
        ///                \"data not available\" or
        ///                \"data not given\".
        /// \remarks If \b S!=NULL at time of call, the function will
        ///  try to dispose the string it points on. This allows a slick
        ///  re-use of the same pointer in consequitive calls. This also
        ///  means that one should never pass unallocated pointer to
        ///  this function. Safe use assumes the following patern:
        ///  \code
        ///  mmcif::Struct mmCIFStruct;
        ///  pstr S;  // this is merely "char *S"
        ///  int  rc;
        ///
        ///    S  = NULL; // null pointer before first use
        ///    rc = mmCIFStruct.GetString ( S,"id" );
        ///    if (rc)  CreateCopy ( S,"*** data not found" );
        ///    if (!S)  CreateCopy ( S,"*** data not given or not available" );
        ///    printf ( " rc=%i, S='%s'\n",rc,S );
        ///
        ///    rc = mmCIFStruct.GetString ( S,"property" );
        ///    if (rc)  CreateCopy ( S,"*** data not found" );
        ///    if (!S)  CreateCopy ( S,"*** data not given or not available" );
        ///    printf ( " rc=%i, S='%s'\n",rc,S );
        ///
        ///    // etc etc etc
        ///
        ///    delete[] S;  // application is responsible for final
        ///                 // disposal of memory
        ///  \endcode
        int  GetString ( pstr & S, cpstr TName, bool Remove=false );

        /// \brief Returns pointer to value associated with given tag.
        /// \param[in] TName character string with tag name
        /// \param[out] RC return code:
        ///    \arg \b CIFRC_NoTag: tag is not found
        ///    \arg \b CIFRC_NoField: value is not found
        ///    \arg \b CIFRC_Ok: success. If function returns NULL, then
        ///                the value corresponds to either
        ///                \"data not available\" or
        ///                \"data not given\".
        /// \return \b NULL: either tag or value is not found, or the
        ///    value is \"data not available\" or \"data not given\".
        ///    Read return code \b RC in order to interpret NULL return.
        /// \return \b not \b NULL: pointer (\c char \c *) to value
        ///    associated with \b TName.
        /// \remarks Never try to dispose memory pointed by the return
        /// value, or your program will crash.
        pstr GetString ( cpstr TName, int & RC ); // NULL if TName
                                                         // is not there

        /// \brief Deletes field associated with given tag.
        /// \param[in] TName character string with tag name
        /// \return \b >=0: field deleted
        /// \return \b <0: either field or tag is not found
        int  DeleteField ( cpstr TName );  // <0 the field was not there

        /// \brief Fetches value, corresponding to the given tag, as
        ///        a real number.
        /// \param[out] R reference to real number to accept the value.
        ///        In case of failure, \b R returns zero.
        /// \param[in] TName character string with tag name
        /// \param[in] Remove flag to remove the tag and its value from
        ///               structure after it is read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_NoField: field is not found
        /// \return \b CIFRC_WrongFormat: value is not a real or integer
        ///            number.
        /// \return \b CIFRC_NoData: value is either
        ///            \"data not available\" or
        ///            \"data not given\".
        /// \return \b CIFRC_Ok: success.
        int  GetReal ( realtype & R, cpstr TName, bool Remove=false );

        /// \brief Fetches value, corresponding to the given tag, as
        ///        an integer number.
        /// \param[out] I reference to integer number to accept the
        ///        value. In case of failure, \b I returns zero, except
        ///        when value is \"data not available\" or
        ///        \"data not given\", when I returns \c MinInt4.
        /// \param[in] TName character string with tag name
        /// \param[in] Remove flag to remove the tag and its value from
        ///               structure after it is read.
        /// \return \arg \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_NoField: field is not found
        /// \return \b CIFRC_WrongFormat: value is not an integer number.
        /// \return \b CIFRC_NoData: value is either
        ///            \"data not available\" or
        ///            \"data not given\".
        /// \return \b CIFRC_Ok: success.
        int  GetInteger ( int & I, cpstr TName, bool Remove=false );

        /// \brief Sets string value for given tag.
        /// \param[in] S character string with value to be set.
        ///            If \b S==NULL, the \"data not given\" value
        ///            will be set. If \b S==\"\" (empty string), the
        ///            \"data not available\" value is stored.
        /// \param[in] TName character string with tag name. If tag
        ///            is not found, it will be added to the structure.
        /// \param[in] NonBlankOnly flag to treat white-space-only
        ///            strings:
        ///   \arg \b false: set as is
        ///   \arg \b true:  set \"data not available\" value instead.
        void PutString   ( cpstr S, cpstr TName,
                           bool NonBlankOnly=false );

        /// \brief Sets current date in format YYYY-MM-DD as a value
        ///        for given tag.
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added to the structure.
        void PutDate     ( cpstr T );

        /// \brief Sets \"data not given\" or \"data not available\"
        ///        values for given tag.
        /// \param[in] NoDataType can be either
        ///   \arg \b CIF_NODATA_DOT for \"data not given\"
        ///   \arg \b CIF_NODATA_QUESTION for \"data not available\"
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added to the structure.
        void PutNoData   ( int NoDataType, cpstr T  );

        /// \brief Sets float-point value for given tag.
        /// \param[in] R real number with value to be set.
        /// \param[in] TName character string with tag name. If tag
        ///            is not found, it will be added to the structure.
        /// \param[in] prec float-point precision; g-format with given
        ///            precision will be used
        void PutReal     ( realtype R, cpstr TName, int prec=8 );

        /// \brief Sets float-point value for given tag.
        /// \param[in] R real number with value to be set.
        /// \param[in] TName character string with tag name. If tag
        ///            is not found, it will be added to the structure.
        /// \param[in] format format string to convert \b R.
        void PutReal     ( realtype R, cpstr TName, cpstr format );

        /// \brief Sets integer value for given tag.
        /// \param[in] I integer number with value to be set.
        /// \param[in] TName character string with tag name. If tag
        ///            is not found, it will be added to the structure.
        void PutInteger  ( int      I, cpstr TName );

        /// \brief Writes structure data in mmCIF format into file.
        /// \param[in] FName character string with file name.
        /// \param[in] gzipMode flag to controll compression of files:
        ///  \arg \b GZM_NONE: do not compress
        ///  \arg \b GZM_CHECK: check file name suffix and compress
        ///                     (or not) accordingly
        ///  \arg \b GZM_ENFORCE_GZIP: force gzip compression despite
        ///                     suffix
        ///  \arg \b GZM_ENFORCE_COMPRESS: force using compress despite
        ///                     suffix
        ///  \return \b true: success
        ///  \return \b false: can not open file for writing.
        /// \remarks This function does not create a valid mmCIF file
        /// as \"data_XXX\" record will be missing. It may be used for
        /// debugging though.
        bool WriteMMCIFStruct ( cpstr FName,
                                io::GZ_MODE gzipMode=io::GZM_CHECK );

        /// \brief Writes structure into given file.
        /// \param f reference to MMDB's file class. The file should be
        /// opened and closed by application.
        /// \remarks There is a very limited use of this function on
        /// application level. It is primarily used by mmcif::Data class.
        void    WriteMMCIF ( io::RFile f  );

        /// \brief Deep copy of structures.
        /// Deep copy duplicates all data and memory allocations,
        /// producing a genuine clone of the original. Only deep copy
        /// should be used for copying MMDB objects, a mere assignment
        /// operator will fail you.
        /// \param[in] Struct a pointer to mmcif::Struct, the content of
        ///                 which is copied into 'this' structure.
        void Copy ( PCategory Struct );

        /// \brief MMDB stream writer.
        void write ( io::RFile f );

        /// \brief MMDB stream reader.
        void read  ( io::RFile f );

      protected:
        psvector field;

        void InitStruct();
        void FreeMemory();

    };



    //  ======================  Loop  ==============================

    DefineClass(Loop);
    DefineStreamFunctions(Loop);

    /// \brief mmcif::Loop represents mmCIF's \"loop\" category, which keeps
    ///        rows of data values associated with tags.
    /*!
    mmCIF's \"loop\" category has the following form:
    \code
    loop_
    _loop_name.tag0   value0
    _loop_name.tag1   value1
    ...........
    _loop_name.tagN   valueN
    value00 value10 ... valueN0
    value01 value11 ... valueN1
    ...........
    value0M value1M ... valueNM
    \endcode
    mmcif::Loop represents this construct by keeping category name
    (\"_loop_name\") and associated lists of tags
    (\"tag0,tag1...tagN\") and data vectors
    (\"[value00...value0M],[value10...value1M]...[valueN0...valueNM]\").

    The loop object is created automatically when an mmCIF file is read,
    or it may be created programatically and then pushed into file.

    Access to data is provided via tags and data indexes. Internally,
    all values are kept as character fields, and it is only on the
    retrieval stage that they are converted to other data types
    (integers, floats or strings). If conversion is not possible, an
    error code is returned by the corresponding functions, which should
    be checked by the application.

    The following code gives an example of creating mmCIF loop category
    and populating it with data:
    \code
    mmcif::Loop loop;
    char       S[100];
    int        i;

      // Specify loop name:
      loop.SetCategoryName ( "_sample_loop" );

      // Create loop structure, i.e., list of tags first:
      loop.AddLoopTag ( "id"    );
      loop.AddLoopTag ( "name"  );
      loop.AddLoopTag ( "value" );

      // Now populate it with data. This my be done in 2 ways.
      // Here is how you write loop data in stream fashion,
      // value-after-value:
      for (i=0;i<3;i++)  {
        loop.AddInteger ( i );
        sprintf ( S,"1st_way-%i",i );
        loop.AddString ( S );
        loop.AddReal ( 2.5*(i+1) );
      }

      // Here is how you populate loop data using direct-access
      // functions:
      for (i=3;i<6;i++)  {
        loop.PutReal ( 2.5*(i+1),"value",i );
        loop.PutInteger ( i,"id" );
        sprintf ( S,"2nd way: %i",i );
        loop.PutString ( S,"name" );
      }

      loop.WriteMMCIFLoop ( "sample_loop.cif" );

    \endcode

    The resulting file \b sample_loop.cif will contain:

    \code

    loop_
    _sample_loop.id
    _sample_loop.name
    _sample_loop.value
    0   1st_way-0     2.5
    1   1st_way-1     5.0
    2   1st_way-2     7.5
    3   "2nd way: 3"  10.0
    4   "2nd way: 4"  12.5
    5   "2nd way: 5"  15.0

    \endcode

    See general principles of working with mmCIF files and mmCIF
    hierarchies, as well as some code samples, in Section
    \"\ref mmcif_handler\".
    */

    class MMDB_IMEX Loop : public Category  {

      friend class Data;

      public :

        /// \brief Basic constructor
        Loop ();

        /// \brief Constructor that assigns structure name.
        /// \param[in] N structure name (must start with underscore).
        Loop ( cpstr N );

        /// \brief Constructor for MMDB data streaming functions
        Loop ( io::RPStream Object );

        /// \brief Destructor
        ~Loop();

        /// \brief Adds tag to the loop.
        /// The tag is appended to the list of existing tags. The order
        /// of tags cannot be changed.
        /// \param[in] T tag name
        /// \param[in] Remove flag to remove all fields in the loop.
        void  AddLoopTag ( cpstr T, bool Remove=true );

        /// \brief Sets string value at current loop position.
        /// When \b mmcif::Loop::Add[Data] functions use internal loop
        /// pointer. When category is created or cleared (by using
        /// \b mmcif::Loop::AddLoopTag ( T,true )) the pointer is set to
        /// 0th row and 0th column (tag). After each call to
        /// \b mmcif::Loop::Add[Data] function, internal pointer advances
        /// to next column (tag), and wraps over to next row, 0th tag,
        /// if list of tags is exhausted. Any remaining fields in last
        /// row will be populated with \"data not given\" value.
        /// \param[in] S character string with value to be set.
        ///            If \b S==NULL, the \"data not given\" value
        ///            will be set. If \b S==\"\" (empty string), the
        ///            \"data not available\" value is stored.
        /// \param[in] NonBlankOnly flag to treat white-space-only
        ///            strings:
        ///   \arg \b false: set as is
        ///   \arg \b true:  set \"data not available\" value instead.
        void  AddString ( cpstr S, bool NonBlankOnly=false );

        /// \brief Sets \"data not given\" or \"data not available\" at
        ///        current loop position.
        /// When \b mmcif::Loop::Add[Data] functions use internal loop
        /// pointer. When category is created or cleared (by using
        /// \b mmcif::Loop::AddLoopTag ( T,true )) the pointer is set to
        /// 0th row and 0th column (tag). After each call to
        /// \b mmcif::Loop::Add[Data] function, internal pointer advances
        /// to next column (tag), and wraps over to next row, 0th tag,
        /// if list of tags is exhausted. Any remaining fields in last
        /// row will be populated with \"data not given\" value.
        /// \param[in] NoDataType integer key specifying which type of
        ///            data absence should be set as a value:
        ///   \arg \b CIF_NODATA_DOT for \"data not given\"
        ///   \arg \b CIF_NODATA_QUESTION for \"data not available\"
        void  AddNoData ( int NoDataType );

        /// \brief Sets float-point value at current loop position.
        /// When \b mmcif::Loop::Add[Data] functions use internal loop
        /// pointer. When category is created or cleared (by using
        /// \b mmcif::Loop::AddLoopTag ( T,true )) the pointer is set to
        /// 0th row and 0th column (tag). After each call to
        /// \b mmcif::Loop::Add[Data] function, internal pointer advances
        /// to next column (tag), and wraps over to next row, 0th tag,
        /// if list of tags is exhausted. Any remaining fields in last
        /// row will be populated with \"data not given\" value.
        /// \param[in] R real number with value to be set.
        /// \param[in] prec float-point precision; g-format with given
        ///            precision will be used
        void  AddReal ( realtype R, int prec=8 );

        /// \brief Sets float-point value at current loop position in
        ///        given format.
        /// When \b mmcif::Loop::Add[Data] functions use internal loop
        /// pointer. When category is created or cleared (by using
        /// \b mmcif::Loop::AddLoopTag ( T,true )) the pointer is set to
        /// 0th row and 0th column (tag). After each call to
        /// \b mmcif::Loop::Add[Data] function, internal pointer advances
        /// to next column (tag), and wraps over to next row, 0th tag,
        /// if list of tags is exhausted. Any remaining fields in last
        /// row will be populated with \"data not given\" value.
        /// \brief Sets float-point value for given tag.
        /// \param[in] R real number with value to be set.
        /// \param[in] format format string to convert \b R.
        void  AddReal ( realtype R, cpstr format );

        /// \brief Sets integer value at current loop position in given
        ///        format.
        /// When \b mmcif::Loop::Add[Data] functions use internal loop
        /// pointer. When category is created or cleared (by using
        /// \b mmcif::Loop::AddLoopTag ( T,true )) the pointer is set to
        /// 0th row and 0th column (tag). After each call to
        /// \b mmcif::Loop::Add[Data] function, internal pointer advances
        /// to next column (tag), and wraps over to next row, 0th tag,
        /// if list of tags is exhausted. Any remaining fields in last
        /// row will be populated with \"data not given\" value.
        /// \param[in] I integer number with value to be set.
        void  AddInteger ( int I );

        /// \brief Returns current length of the loop (i.e. the number
        ///        of rows).
        /// \return number of data rows in the loop.
        int   GetLoopLength() { return nRows; }

        /// \brief Returns string pointer on the field corresponding to
        ///        tag in the specified position, in the specified row.
        /// Tag positions are defined by the order of their appearance in
        /// mmCIF file (if loop was read from a file), or by the
        /// order of their addition to the loop (if loop was
        /// created programmatically).
        /// \param[in] rowNo row number (0...GetLoopLength()-1)
        /// \param[in] tagNo tag number (0...GetNofTags()-1)
        /// \return \b NULL: tag or row do not exist
        /// \return \b CIF_NODATA_DOT_FIELD the field contains
        ///            \"data not given\" value
        /// \return \b CIF_NODATA_QUESTION_FIELD the field contains
        ///            \"data not available\" value
        /// \return \b not \b NULL: string value of the field
        /// \remarks Never try to dispose memory pointed by the return
        /// value, or your program will crash.
        pstr  GetField ( int rowNo, int tagNo );

        /// \brief Fetches value, corresponding to the given tag, in
        ///        the given row, as a string
        /// \param[out] S pointer to string, which will point to newly
        ///               allocated character string, containing value
        ///               associated with tag \b TName and row \b nrow.
        ///               If tag, row or value
        ///               is not found, or if value corresponds to
        ///               mmCIF's \"data not given\" or
        ///               \"data not available\", \b S returns NULL.
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[in] Remove flag to remove the field from
        ///               structure after it is read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_WrongIndex: row is not found
        /// \return \b CIFRC_NoField: value is not found
        /// \return \b CIFRC_Ok: success. If \b S returns NULL, then
        ///                the value corresponds to either
        ///                \"data not available\" or
        ///                \"data not given\".
        /// \remarks If \b S!=NULL at time of call, the function will
        ///  try to dispose the string it points on. This allows a slick
        ///  re-use of the same pointer in consequitive calls. This also
        ///  means that one should never pass unallocated pointer to
        ///  this function. Safe use assumes the following patern:
        ///  \code
        ///  mmcif::Loop mmCIFLoop;
        ///  pstr S;  // this is merely "char *S"
        ///  int  rc;
        ///
        ///    S  = NULL; // null pointer before first use
        ///    rc = mmCIFLoop.GetString ( S,"id",1 );
        ///    if (rc)  CreateCopy ( S,"*** data not found" );
        ///    if (!S)  CreateCopy ( S,"*** data not given or not available" );
        ///    printf ( " rc=%i, S='%s'\n",rc,S );
        ///
        ///    rc = mmCIFLoop.GetString ( S,"property",0 );
        ///    if (rc)  CreateCopy ( S,"*** data not found" );
        ///    if (!S)  CreateCopy ( S,"*** data not given or not available" );
        ///    printf ( " rc=%i, S='%s'\n",rc,S );
        ///
        ///    // etc etc etc
        ///
        ///    delete[] S;  // application is responsible for final
        ///                 // disposal of memory
        ///  \endcode
        int   GetString ( pstr & S, cpstr TName, int nrow,
                                           bool Remove=false );

        /// \brief Returns pointer to value associated with given tag,
        ///        in the given row of the loop.
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[out] RC return code:
        ///    \arg \b CIFRC_NoTag: tag is not found
        ///    \arg \b CIFRC_WrongIndex: row is not found
        ///    \arg \b CIFRC_NoField: value is not found
        ///    \arg \b CIFRC_Ok: success. If function returns NULL, then
        ///                the value corresponds to either
        ///                \"data not available\" or
        ///                \"data not given\".
        /// \return \b NULL: either tag, row or value is not found, or the
        ///    value is \"data not available\" or \"data not given\".
        ///    Read return code \b RC in order to interpret NULL return.
        /// \return \b not \b NULL: pointer (\c char \c *) to value
        ///    associated with \b TName.
        /// \remarks Never try to dispose memory pointed by the return
        /// value, or your program will crash.
        pstr  GetString    ( cpstr TName, int nrow, int & RC );

        /// \brief Copies value, associated with given tag,
        ///        in the given row, into specified buffer.
        ///  Terminating NULL character is appended.
        /// \param[out] buf character string to accept the value
        /// \param[in] maxlength maximum number of bytes to copy
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[out] RC return code:
        ///    \arg \b CIFRC_NoTag: tag is not found
        ///    \arg \b CIFRC_WrongIndex: row is not found
        ///    \arg \b CIFRC_NoField: value is not found
        ///    \arg \b CIFRC_Ok: success.
        /// \remarks Destination string \b buf is not modified if
        /// \b RC!=CIFRC_Ok .
        void  CopyString   ( pstr  buf,   int maxlength,
                             cpstr TName, int nrow, int & RC );

        /// \brief Deletes field associated with given tag in
        ///          the given row.
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \return \b >=0: field deleted
        /// \return \b <0: either field or tag is not found
        int   DeleteField  ( cpstr TName, int nrow );

        /// \brief Deletes all fields in given row.
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \return \b CIFRC_Ok: fields deleted
        /// \return \b CIFRC_WrongIndex: row not found
        /// \remarks Note that this function delets just the fields, but
        /// not the row. If you wish the row to be deleted, call
        /// mmcif::Loop::Optimize() function after this one.
        int   DeleteRow    ( int nrow );

        /// \brief Fetches value, corresponding to the given tag,
        ///        in the given row, as a real number.
        /// \param[out] R reference to real number to accept the value.
        ///        In case of failure, \b R returns zero.
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[in] Remove flag to remove the field from
        ///               the loop after it is read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_WrongIndex: row not found
        /// \return \b CIFRC_NoField: field is not found
        /// \return \b CIFRC_WrongFormat: value is not a real or integer
        ///            number.
        /// \return \b CIFRC_NoData: value is either
        ///            \"data not available\" or
        ///            \"data not given\".
        /// \return \b CIFRC_Ok: success.
        int   GetReal ( realtype & R, cpstr TName, int nrow,
                                           bool Remove=false );

        /// \brief Copies value, associated with given tag,
        ///        in the given row, into specified destination as
        ///        a real number.
        /// \param[out] R reference to real number to accept the value
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[out] RC return code:
        ///    \arg \b CIFRC_NoTag: tag is not found
        ///    \arg \b CIFRC_WrongIndex: row is not found
        ///    \arg \b CIFRC_NoField: value is not found
        ///    \arg \b CIFRC_Ok: success.
        /// \remarks Destination \b R is set 0 if \b RC!=CIFRC_Ok.
        void  CopyReal ( realtype & R, cpstr TName, int nrow, int & RC );

        /// \brief Copies value, associated with given tag,
        ///        in the given row, into specified destination as
        ///        an integer number.
        /// \param[out] I reference to integer number to accept the value
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[out] RC return code:
        ///    \arg \b CIFRC_NoTag: tag is not found
        ///    \arg \b CIFRC_WrongIndex: row is not found
        ///    \arg \b CIFRC_NoField: value is not found
        ///    \arg \b CIFRC_Ok: success.
        /// \remarks Destination \b I is set 0 if \b RC!=CIFRC_Ok.
        void  CopyInteger ( int & I, cpstr TName, int nrow, int & RC );

        /// \brief Fetches value, corresponding to the given tag,
        ///        in the given row, as an integer number.
        /// \param[out] I reference to integer number to accept the value.
        ///        In case of failure, \b R returns zero.
        /// \param[in] TName character string with tag name
        /// \param[in] nrow  row number (0...GetLoopLength()-1)
        /// \param[in] Remove flag to remove the field from
        ///               the loop after it is read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_WrongIndex: row not found
        /// \return \b CIFRC_NoField: field is not found
        /// \return \b CIFRC_WrongFormat: value is not a real or integer
        ///            number.
        /// \return \b CIFRC_NoData: value is either
        ///            \"data not available\" or
        ///            \"data not given\".
        /// \return \b CIFRC_Ok: success.
        int   GetInteger   ( int & I, cpstr TName, int nrow,
                                           bool Remove=false );

        /// \brief Fetches set of values, corresponding to the given
        ///        tag, in the given range of rows, as a vector of
        ///        strings.
        /// \param[out] S reference to string vector to accept
        ///        the values. if \b S==NULL , the vector will be
        ///        allocated with starting index of \b i1.
        /// \param[in] TName character string with tag name
        /// \param[in] i1  minimum row number to fetch, the actual
        ///            index will be calculated as \b max(0,min(i1,i2))
        /// \param[in] i2  maximum row number to fetch, the actual
        ///            index will be calculated as
        ///            \b min(GetLoopLength()-1,max(i1,i2))
        /// \param[in] Remove flag to remove fetched fields from
        ///               the loop after they are read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_WrongIndex: invalid range of rows
        /// \return \b CIFRC_Ok: success.
        ///
        /// For safe use, \b S should be pre-allocated by calling
        /// process. Only elements \b S[i1] to \b S[i2] will contain
        /// fetched data, others remain untouched. The calling
        /// process is responsible for the disposal of \b S. Example:
        /// \code
        /// mmcif::Loop loop;
        /// psvector   S;  // equivalent to char **S
        /// int        i,i1,i2,rc,n;
        ///
        ///    // ... get loop data
        ///
        ///    n  = loop.GetLoopLength();
        ///    i1 = 5;  i2 = n - 5;  // could be wrong!
        ///
        ///    //  allocate vector of strings
        ///    GetVectorMemory ( S,n,0 );  // "0" for starting index
        ///    for (i=0;i<n;i++)
        ///      S[i] = NULL;  // initialize NULL string pointers
        ///
        ///    loop.GetSVector ( S,"name",i1,i2 );
        ///    printf ( " Fetched with return code rc=%i\n",rc );
        ///        // you may want a more thorough treatment of
        ///        // the return code here
        ///    for (i=i1;i<=i2;i++)
        ///      if (S[i])  printf ( " %4i. name='%s'\n",i,S[i] );
        ///           else  printf ( " %4i. name is not available\n",i );
        ///
        ///    //  S[] may be re-used for as many fetches as necessary
        ///    //  without cleaning or disposals
        ///
        ///    //  dispose of vector of strings
        ///    for (i=0;i<n;i++)
        ///      if (S[i])  delete[] S[i];
        ///    FreeVectorMemory ( S,0 );  // "0" for starting index
        ///
        /// \endcode
        int  GetSVector ( psvector & S, cpstr TName,
                          int i1=0, int i2=MaxInt4,
                          bool Remove=false );

        /// \brief Fetches set of values, corresponding to the given
        ///        tag, in the given range of rows, as a vector of
        ///        float-point numbers.
        /// \param[out] R reference to float-point vector to accept
        ///        the values. if \b R==NULL , the vector will be
        ///        allocated with starting index of \b i1.
        /// \param[in] TName character string with tag name
        /// \param[in] i1  minimum row number to fetch, the actual
        ///            index will be calculated as \b max(0,min(i1,i2))
        /// \param[in] i2  maximum row number to fetch, the actual
        ///            index will be calculated as
        ///            \b min(GetLoopLength()-1,max(i1,i2))
        /// \param[in] Remove flag to remove fetched fields from
        ///               the loop after they are read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_WrongIndex: invalid range of rows
        /// \return \b CIFRC_Ok: success.
        ///
        /// For safe use, \b R should be pre-allocated by calling
        /// process. Only elements \b R[i1] to \b R[i2] will contain
        /// fetched data, others remain untouched. The calling
        /// process is responsible for the disposal of \b R. Example:
        /// \code
        /// mmcif::Loop loop;
        /// rvector    R;  // equivalent to realtype *R
        /// int        i,i1,i2,rc,n;
        ///
        ///    // ... get loop data
        ///
        ///    n  = loop.GetLoopLength();
        ///    i1 = 5;  i2 = n - 5;  // could be wrong!
        ///
        ///    //  allocate a vector of real numbers
        ///    GetVectorMemory ( R,n,0 );  // "0" for starting index
        ///    // no need to initiaize unless required for the
        ///    // application
        ///
        ///    rc = loop.GetRVector ( R,"value",i1,i2 );
        ///    printf ( " Fetched with return code rc=%i\n",rc );
        ///        // you may want a more thorough treatment of
        ///        // the return code here
        ///    for (i=i1;i<=i2;i++)
        ///      printf ( " value[%4i] = %15.7g\n",i,R[i] );
        ///
        ///    //  R[] may be re-used for as many fetches as necessary
        ///    //  without cleaning or disposals
        ///
        ///    //  dispose of the vector
        ///    FreeVectorMemory ( R,0 );  // "0" for starting index
        ///
        /// \endcode
        int  GetRVector ( rvector  & R, cpstr TName,
                          int i1=0, int i2=MaxInt4,
                          bool Remove=false );

        /// \brief Fetches set of values, corresponding to the given
        ///        tag, in the given range of rows, as a vector of
        ///        integer numbers.
        /// \param[out] I reference to float-point vector to accept
        ///        the values. if \b I==NULL , the vector will be
        ///        allocated with starting index of \b i1.
        /// \param[in] TName character string with tag name
        /// \param[in] i1  minimum row number to fetch, the actual
        ///            index will be calculated as \b max(0,min(i1,i2))
        /// \param[in] i2  maximum row number to fetch, the actual
        ///            index will be calculated as
        ///            \b min(GetLoopLength()-1,max(i1,i2))
        /// \param[in] Remove flag to remove fetched fields from
        ///               the loop after they are read.
        /// \return \b CIFRC_NoTag: tag is not found
        /// \return \b CIFRC_WrongIndex: invalid range of rows
        /// \return \b CIFRC_Ok: success.
        ///
        /// For safe use, \b I should be pre-allocated by calling
        /// process. Only elements \b I[i1] to \b I[i2] will contain
        /// fetched data, others remain untouched. The calling
        /// process is responsible for the disposal of \b I.
        /// See example in mmcif::Loop::GetRVector documentation
        /// for details.
        int  GetIVector ( ivector  & I, cpstr TName,
                          int i1=0, int i2=MaxInt4,
                          bool Remove=false );

        /// \brief Sets string value for given tag and row.
        /// \param[in] S character string with value to be set.
        ///            If \b S==NULL, the \"data not given\" value
        ///            will be set. If \b S==\"\" (empty string), the
        ///            \"data not available\" value is stored.
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] nrow  row number. If the row does not exist then
        ///            it will be created, along with all other rows
        ///            between GetLoopLength()-1 and \b nrow as
        ///            necessary. All newly created fields will be
        ///            initialised with \"data not given\" value.
        void  PutString ( cpstr S, cpstr T, int nrow );

        /// \brief Sets \"data not given\" or \"data not available\"
        ///        values for given tag and row.
        /// \param[in] NoDataType can be either
        ///   \arg \b CIF_NODATA_DOT for \"data not given\"
        ///   \arg \b CIF_NODATA_QUESTION for \"data not available\"
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] nrow  row number. If the row does not exist then
        ///            it will be created, along with all other rows
        ///            between GetLoopLength()-1 and \b nrow as
        ///            necessary. All newly created fields will be
        ///            initialised with \"data not given\" value.
        void  PutNoData ( int NoDataType, cpstr T, int nrow );

        /// \brief Sets float-point value for given tag and row.
        /// \param[in] R real number with value to be set.
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] nrow  row number. If the row does not exist then
        ///            it will be created, along with all other rows
        ///            between GetLoopLength()-1 and \b nrow as
        ///            necessary. All newly created fields will be
        ///            initialised with \"data not given\" value.
        /// \param[in] prec float-point precision; g-format with given
        ///            precision will be used
        void  PutReal ( realtype R, cpstr T, int nrow, int prec=8 );

        /// \brief Sets float-point value for given tag and row.
        /// \param[in] R real number with value to be set.
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] nrow  row number. If the row does not exist then
        ///            it will be created, along with all other rows
        ///            between GetLoopLength()-1 and \b nrow as
        ///            necessary. All newly created fields will be
        ///            initialised with \"data not given\" value.
        /// \param[in] format format string to convert \b R.
        void  PutReal ( realtype R, cpstr T, int nrow, cpstr format );

        /// \brief Sets integer value for given tag.
        /// \param[in] I integer number with value to be set.
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] nrow  row number. If the row does not exist then
        ///            it will be created, along with all other rows
        ///            between GetLoopLength()-1 and \b nrow as
        ///            necessary. All newly created fields will be
        ///            initialised with \"data not given\" value.
        void  PutInteger ( int I, cpstr T, int nrow );

        /// \brief Sets a set of string values for the given tag and
        ///        range of rows.
        /// \param[in] S string vector with values to store in the loop
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] i1  minimum data index in \b S to set in the loop
        /// \param[in] i2  maximum data index in \b S to set in the loop.
        ///
        /// The data will be set in rows \b i1 to \b i2 (inclusive) in
        /// the loop. If range \b [i1,i2] is not contained in the loop,
        /// all missing rows will be created and initialised to
        /// \"data not given\" value. Example:
        /// \code
        /// mmcif::Loop loop("_sample_loop");
        /// pstr       S[100];
        /// int        i;
        ///
        ///    //  initialize vector of strings
        ///    for (i=0;i<100;i++)  {
        ///      S[i] = new char[20];
        ///      sprintf ( S[i],"value i=%i",i );
        ///    }
        ///
        ///    //  put data in loop
        ///    loop.PutSVector ( S,"made_up_string_value",0,99 );
        ///
        ///    //  dispose of vector of strings
        ///    for (i=0;i<100;i++)
        ///      if (S[i])  delete[] S[i];
        ///
        /// \endcode
        void  PutSVector   ( psvector S, cpstr T, int i1, int i2 );

        /// \brief Sets a set of float-point values for the given tag and
        ///        range of rows.
        /// \param[in] R vector of real numbers to store in the loop
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] i1  minimum data index in \b S to set in the loop
        /// \param[in] i2  maximum data index in \b S to set in the loop
        /// \param[in] prec float-point precision; g-format with given
        ///            precision will be used.
        ///
        /// The data will be set in rows \b i1 to \b i2 (inclusive) in
        /// the loop. If range \b [i1,i2] is not contained in the loop,
        /// all missing rows will be created and initialised to
        /// \"data not given\" value.
        void  PutRVector   ( rvector  R, cpstr T, int i1, int i2,
                                                       int prec=8 );

        /// \brief Sets a set of integer values for the given tag and
        ///        range of rows.
        /// \param[in] I vector of integers to store in the loop
        /// \param[in] T character string with tag name. If tag
        ///            is not found, it will be added, and all data in
        ///            the loop will be reindexed accordingly.
        /// \param[in] i1  minimum data index in \b S to set in the loop
        /// \param[in] i2  maximum data index in \b S to set in the loop.
        ///
        /// The data will be set in rows \b i1 to \b i2 (inclusive) in
        /// the loop. If range \b [i1,i2] is not contained in the loop,
        /// all missing rows will be created and initialised to
        /// \"data not given\" value.
        void  PutIVector   ( ivector  I, cpstr T, int i1, int i2 );

        /// \brief Returns category type \b MMCIF_Loop.
        MMCIF_ITEM  GetCategoryID() { return MMCIF_Loop; }

        /// \brief Optimizes loop for RAM and data access speed.
        /// Optimized data structures take less RAM and their indexes
        /// are sorted for quicker access. Sorting is done automatically
        /// as new data is added to the category. If the structure
        /// is edited (fields/data removed), it may need
        /// optimization and re-sorting for efficiency.\n\n
        /// The sorting preserves the order of actual appearance of
        /// tags and rows in mmCIF file. If a loop is created
        /// programmatically, the order of tags and rows in mmCIF file
        /// will be the same as order of adding them to the loop.
        void  Optimize();

        /// \brief Writes loop data in mmCIF format into file.
        /// \param[in] FName character string with file name.
        /// \param[in] gzipMode flag to controll compression of files:
        ///  \arg \b GZM_NONE: do not compress
        ///  \arg \b GZM_CHECK: check file name suffix and compress
        ///                     (or not) accordingly
        ///  \arg \b GZM_ENFORCE_GZIP: force gzip compression despite
        ///                     suffix
        ///  \arg \b GZM_ENFORCE_COMPRESS: force using compress despite
        ///                     suffix
        /// \return \b true: success
        /// \return \b false: can not open file for writing.
        /// \remarks This function does not create a valid mmCIF file
        /// as \"data_XXX\" record will be missing. It may be used for
        /// debugging though.
        bool WriteMMCIFLoop ( cpstr FName,
                              io::GZ_MODE gzipMode=io::GZM_CHECK );

        /// \brief Writes loop data into given file.
        /// \param f reference to MMDB's file class. The file should be
        /// opened and closed by application.
        /// \remarks There is a very limited use of this function on
        /// application level. It is primarily used by mmcif::Data class.
        void  WriteMMCIF ( io::RFile f );

        /// \brief Deep copy of loops.
        /// Deep copy duplicates all data and memory allocations,
        /// producing a genuine clone of the original. Only deep copy
        /// should be used for copying MMDB objects, a mere assignment
        /// operator will fail you.
        /// \param[in] Loop a pointer to mmcif::Loop, the content of
        ///                 which is copied into 'this' loop.
        void  Copy ( PCategory Loop );

        /// \brief MMDB stream writer.
        void write ( io::RFile f );

        /// \brief MMDB stream reader.
        void read  ( io::RFile f );

      protected:
        int      nRows;
        psmatrix field;
        int      iColumn,nAllocRows;

        void  InitLoop     ();
        void  FreeMemory   ();
        void  DeleteFields ();
        void  ExpandRows   ( int nRowsNew );

    };



    //  ======================  Data  =============================


    //    CIFW are warnings which may be issued on reading the CIF file.
    // Each of them means actually a CIF syntax error.

    enum CIF_WARNING  {
      CIFW_UnrecognizedItems = 0x00000020,
      CIFW_MissingField      = 0x00000040,
      CIFW_EmptyLoop         = 0x00000080,
      CIFW_UnexpectedEOF     = 0x00000100,
      CIFW_LoopFieldMissing  = 0x00000200,
      CIFW_NotAStructure     = 0x00000400,
      CIFW_NotALoop          = 0x00000800,
      CIFW_DuplicateTag      = 0x00001000
    };

    //    CIFRC are return codes from procedures of extracting data from
    // the read CIF file. Negative returns reflect unsuccessful and
    // not accomplished operation.
    enum CIF_RC  {
      CIFRC_Loop           =  2,
      CIFRC_Structure      =  1,
      CIFRC_Ok             =  0,
      CIFRC_StructureNoTag = -1,
      CIFRC_LoopNoTag      = -2,
      CIFRC_NoCategory     = -3,
      CIFRC_WrongFormat    = -4,
      CIFRC_NoTag          = -5,
      CIFRC_NotAStructure  = -6,
      CIFRC_NotALoop       = -7,
      CIFRC_WrongIndex     = -8,
      CIFRC_NoField        = -9,
      CIFRC_Created        = -12,
      CIFRC_CantOpenFile   = -13,
      CIFRC_NoDataLine     = -14,
      CIFRC_NoData         = -15
    };

    //
    //  Functional flags:
    //  ~~~~~~~~~~~~~~~~~
    //
    //  CIFFL_PrintWarnings      when reading CIF file, all warning
    //                           messages will be printed. If the flag
    //                           is off, the warnings will be bit-encoded
    //                           in the return code
    //  CIFFL_StopOnWarnings     reading CIF file will stop at first
    //                           warning issued
    //  CIFFL_SuggestCategories  allows reading CIF file with loops having
    //                           no categories. Hidden category names
    //                           will be automatically generated for
    //                           internal consistency of the system.
    //                           These names will not appear in output.
    //                           As these names are hidden, they cannot
    //                           be used to access data. It is therefore
    //                           assumed that all tags in all loops without
    //                           categories are unique. Simply specify ""
    //                           for category when accessing such data
    //                           (it cannot be accessed through mmcif::Loop,
    //                           but only through mmcif::Data functions
    //                           taking both Category and Tag; note that
    //                           CIFFL_SuggestCategories flag must be on
    //                           while accessing such data).
    //  CIFFL_SuggestTags        allows for identical tags in a category
    //                           (including a hidden category). Hidden
    //                           suffixes to tag names will be generated
    //                           for internal consistency. At present,
    //                           only data for first non-unique tag may be
    //                           accessed.
    //
    enum CIF_FLAG  {
      CIFFL_PrintWarnings     = 0x00000001,
      CIFFL_StopOnWarnings    = 0x00000002,
      CIFFL_SuggestCategories = 0x00000004,
      CIFFL_SuggestTags       = 0x00000008
    };

    DefineClass(Data);
    DefineStreamFunctions(Data);


    /// \brief mmcif::Data represents mmCIF's \"data\" category, which keeps
    ///        structures and loops and is mandatory element of mmCIF file.
    /*!
    mmCIF's \"data\" category has the following form:
    \code
    data_DataName

    _structure1.tag1  value1
    ..........

    loop_
    ..........

    \endcode
    In the above example, all structures and loops that follow \b data_
    keyword until next \b data_ or end of file are part of data category
    with name \b DataName.

    mmcif::Data represents this construct by keeping a list of mmcif::Struct
    and mmcif::Loop class instances associated with the corresponding
    categories in the data block.

    The data object is created automatically when an mmCIF file is read,
    or it may be created programatically and then pushed into file.

    Access to data is provided via category (structures and loops) names,
    tags and data indexes (in case of loops). Alternatively, pointers to
    contained structures and loops may be obtained first, an used for
    fetching data using mmcif::Struct's and mmcif::Loop's interface
    functions.

    The following code gives an example of creating mmCIF's data category
    and populating it:
    \code
    mmcif::Data data;

      // Specify data name:
      data.PutDataName ( "Sample_Data" );

      // the following statement:
      data.PutInteger ( 12345,"_category1","id" );
      // creates structure "_category1" with tag "id" and assigns it
      // the integer value of 12345.

      data.PutString ( "a name","_category1","name" );

      // Loops may be created quite similarly:
      data.PutLoopInteger ( 12345   ,"_loop1","id"  ,2 );
      data.PutLoopInteger ( "a name","_loop1","name",0 );

      // push data into a file
      data.WriteMMCIFData ( "sample.cif" );

    \endcode

    The resulting file \b sample.cif will contain:

    \code
    data_Sample_Data

    _category1.id   12345
    _category1.name "a name"

    loop_
    _loop1.id
    _loop1.name
    .      "a name"
    .      .
    12345  .
    \endcode

    The same result may be achieved differently:

    \code
    mmcif::Data    data;
    mmcif::PStruct mmCIFStruct;  // equivalent to mmcif::Struct *mmCIFStruct
    mmcif::PLoop   mmCIFLoop;    // equivalent to mmcif::Loop   *mmCIFLoop

      // Specify data name:
      data.PutDataName ( "Sample_Data" );

      // create new mmCIF's structure in the data block:
      data.AddStructure ( "_category1",mmCIFStruct );
      if (mmCIFStruct)  {
        mmCIFStruct->PutInteger ( 12345   ,"id"   );
        mmCIFStruct->PutString  ( "a name","name" );
      }

      // similarly for the loop:
      data.AddLoop ( "_loop1",mmCIFLoop );
      if (mmCIFLoop)  {
        mmCIFLoop->PutInteger ( 12345   ,"id"  ,2 );
        mmCIFLoop->PutString  ( "a name","name",0 );
      }

      // push data into a file
      data.WriteMMCIFData ( "sample.cif" );

    \endcode

    See general principles of working with mmCIF files and mmCIF
    hierarchies, as well as some code samples, in Section
    \"\ref mmcif_handler\".
    */

    class MMDB_IMEX Data : public io::Stream  {

      friend class File;

      public :

        /// \brief Basic constructor.
        Data ();

        /// \brief Constructor that assigns data block name.
        /// \param[in] N data block name.
        Data ( cpstr N );

        /// \brief Constructor for MMDB data streaming functions.
        Data ( io::RPStream Object );

        /// \brief Destructor.
        ~Data();


        // -------- General I/O functions

        /// \brief Sets flag to print warnings on reading mmCIF files.
        /// \param[in] SPW flag to print warnings:
        ///    \arg \b true : warnings will be printed to stdout
        ///    \arg \b false : warnings will not be printed but returned
        ///                    in return code (default)
        void  SetPrintWarnings ( bool SPW );

        /// \brief Sets flag to stop on warning when reading an mmCIF file.
        /// \param[in] SOW flag to stop on warning:
        ///    \arg \b true : reading will stop on first warning encountered
        ///    \arg \b false : warnings will not stop reading (default)
        void  SetStopOnWarning ( bool SOW );

        /// \brief Sets optional flag(s) for reading mmCIF files.
        /// By default, no flags are set.
        /// \param[in] F flag or logical \"or\" of several flags to be set:
        ///  \arg \b CIFFL_PrintWarnings  toggles printing warning messages
        ///               at reading an mmCIF file, in stdout. If this
        ///               flag is not set (default), the warnings will
        ///               be returned in the bit-encoded return code
        ///  \arg \b CIFFL_StopOnWarnings  if set, reading an mmCIF file
        ///               will stop at first warning issued
        ///  \arg \b CIFFL_SuggestCategories  allows for reading of mmCIF
        ///               files with loops and structures having no
        ///               category names (\"dirty CIFs\"). If this flag is
        ///               set, then hidden category names will be
        ///               automatically generated. These names will not
        ///               appear in the output. As these names are hidden,
        ///               they cannot be used to access data. In order to
        ///               access data in such categories, consider whether
        ///               they are structures or loops. In case of a
        ///               unnamed structure, simply specify \"\" (empty
        ///               string) for structure name in all access
        ///               functions ( note that \b CIFFL_SuggestCategories
        ///               flag must be on while accessing such data). In
        ///               case of a loop, first use the mmcif::Data::FindLoop
        ///               function to retrieve pointer on the hidden loop,
        ///               and then use mmcif::Loop's interface function to
        ///               fetch data from the loop.
        ///  \arg \b CIFFL_SuggestTags  allows for duplicate tags in a
        ///               category (structure or loop, including hidden
        ///               categories). This may help reading \"dirty CIFs\".
        ///               At present, only data for first non-unique tag
        ///               may be accessed.
        void  SetFlag ( CIF_FLAG F );

        /// \brief Removes optional flag(s) for reading mmCIF files.
        /// By default, no flags are set.
        /// \param[in] F flag or logical \"or\" of several flags to be
        ///              removed (unset):
        ///  \arg \b CIFFL_PrintWarnings  no wornings will be printed in
        ///               stdout, but rather returned in the bit-encoded
        ///               return code
        ///  \arg \b CIFFL_StopOnWarnings  warnings will not stop reading
        ///  \arg \b CIFFL_SuggestCategories  loops without names will
        ///               count as errors and stop reading
        ///  \arg \b CIFFL_SuggestTags  duplicate tags in structures and
        ///               loops will count as errors and stop reading.
        ///
        /// See more detail flag description in mmcif::Data::SetFlag().
        void  RemoveFlag ( CIF_FLAG F );

        /// \brief Returns bit-encoded warnings issued at last file read.
        /// \return an integer number, which is an or-superposition of
        ///         warning flags:
        /// \arg \b CIFW_UnrecognizedItems: unrecognized items were found
        /// \arg \b CIFW_MissingField: expected data field not found
        /// \arg \b CIFW_EmptyLoop: loop category was defined but has no
        ///                         data
        /// \arg \b CIFW_UnexpectedEOF: mmCIF construct finished prematurely
        /// \arg \b CIFW_LoopFieldMissing: loop category has wrong number
        ///                         of data fields
        /// \arg \b CIFW_NotAStructure: attempt to use a category name,
        ///                         which was once defined as a structure,
        ///                         as a loop
        /// \arg \b CIFW_NotALoop: attempt to use a category name, which was
        ///                         once defined as a loop, as a structure
        /// \arg \b CIFW_DuplicateTag: duplicate tag was found in a
        ///                         structure or loop
        inline int  GetWarnings() { return Warning; }

        /// \brief Sets category names and tags that are to be ignored
        ///        on file read.
        /// \param[in] cats list of categories, terminated by NULL
        /// \param[in] tags list of tags, terminated by NULL.
        ///
        /// This special function is to aid reading corrupt mmCIF files.
        /// The input lists should be of equal length 'n', and specify
        /// 'n' \"wrong fields\" that should be ignored on input. E.g.,
        /// ith \"wrong field\" is identified as \"cats[i].taga[i]\".
        /// If \"wrong field\" belongs to a loop, then all the corresponding
        /// column is assumed to be absent. This corrects for mmCIF errors
        /// when defined tags in loops or structures do not have actual data
        /// associated with them.
        ///
        /// In order to remove settings, call SetWrongFields(NULL,NULL).
        ///
        /// Example:
        /*!
        \code
        // assume data for "_category.item1" and "_category.item2"
        // missed in a file to-be-read
        mmcif::Data data;
        cpstr cats[] = { "_category", "_category", NULL };
        cpstr tags[] = { "item1"    , "item2"    , NULL };

           data.SetWrongFields ( cats,tags );
           data.ReadMMCIFData  ( "corrupt.cif" );
        \endcode
        */
        void  SetWrongFields ( cpstr *cats, cpstr *tags );

        /// \brief Reads mmCIF data block from file.
        /// \param FName character null-terminated string with file name
        /// \param gzipMode flag to read compressed files:
        /// \arg \b GZM_NONE: do not assume any compression
        /// \arg \b GZM_CHECK: check compression type by file extension
        /// \arg \b GZM_ENFORCE: same as \b GZM_ENFORCE_GZIP
        /// \arg \b GZM_ENFORCE_GZIP: assume gzip compression (*.gz files)
        /// \arg \b GZM_ENFORCE_COMPRESS: assume compression with 'compress'
        ///         (*.Z files).
        /// \return \b CIFRC_Ok: no errors
        /// \return \b negative: there were errors
        /// \return \b positive: there were warnings.
        ///
        /// This function will read 1st data block from the specified file.
        /// In case of non-zero return, use GetCIFMessage() function to
        /// print the corresponding error message or warning:
        /*!
        \code
        mmcif::Data data;
        char       errLog[500];
        int        rc;
           rc = data.ReadMMCIFData  ( "myfile.cif" );
           if (rc<0)
             printf ( " There was an error:\n %s\n",
                      GetCIFMessage(errLog,rc) );
           else if (rc>0)
             printf ( " There were warnings:\n %s\n",
                      GetCIFMessage(errLog,rc) );
           else
             printf ( " mmCIF file has be read in successfully.\n" );
        \endcode
        */
        int  ReadMMCIFData ( cpstr FName,
                             io::GZ_MODE gzipMode=io::GZM_CHECK );

        /// \brief Reads sequential mmCIF data blocks from file.
        /// \param RCFile reference to a CFile object opened on a file
        /// \param S buffer string which represent a sliding read window.
        ///          The string should be at least 500 characters long,
        ///          initialized with empty-string value before first read,
        ///          and passed unchanged between the reads
        /// \param lcount line counter, should be set zero before first
        ///          read and passed unchanged between the reads.
        /// \return \b CIFRC_Ok: no errors
        /// \return \b negative: there were errors
        /// \return \b positive: there were warnings.
        ///
        /// This function will read 1st data block from the current position
        /// of the file. The function is useful if a file contains more than
        /// a single data block, which should be read sequentially.
        ///
        /// \note Alternatively, files with multiple data blocks can be
        /// read using mmcif::File class.
        ///
        /// In case of non-zero return, use GetCIFMessage() function to
        /// print the corresponding error message or warning:
        /*!
        \code
      mmcif::Data    mmCIFData;
      CFile         f;
      char          S[1000];
      int           rc,lcount;

        // open file first
        f.assign ( "/path/example.cif" );
        if (!f.reset(true))  {
          printf ( " *** cannot open file '%s' for reading.\n",
                   f.FileName() );
          return -1;
        }

        lcount = 0;         // global line counter through the file
        S[0]   = char(0);   // buffer string
        while (!f.FileEnd())  {

          rc = mmCIFData.ReadMMCIFData ( f,S,lcount );

          if (rc!=CIFRC_Ok)  {  // error or warning
            if ((rc<0) && (!f.FileEnd()))  { // error
              printf ( " *** error reading file %s:\n"
                       "     %s\n",f.FileName(),GetCIFMessage(S,rc) );
              return rc;
            } else if (rc>0)  { // warning
              printf ( " ... warning on reading file %s:\n"
                       "     %s\n",f.FileName(),GetCIFMessage(S,rc) );
            }
          } else  {
            // fetch needful values from the data block
            // ........
          }

        }

        f.shut();  // close file

        // NOTE: do not delete mmcif::Struct/mmcif::Loop
        // classes obtained from mmcif::Data. If you do, get a crash.
        // All these structures are containers that dispose their
        // content automatically.
        \endcode
        */
        int  ReadMMCIFData ( io::RFile f, pstr S, int & lcount );

        /// \brief Writes mmCIF data block into file.
        /// \param FName character null-terminated string with file name
        /// \param gzipMode flag to read compressed files:
        /// \arg \b GZM_NONE: do not compress
        /// \arg \b GZM_CHECK: compress according to file extension
        /// \arg \b GZM_ENFORCE: same as \b GZM_ENFORCE_GZIP
        /// \arg \b GZM_ENFORCE_GZIP: compress with gzip
        /// \arg \b GZM_ENFORCE_COMPRESS: compression with 'compress'.
        /// \return \b true: no errors
        /// \return \b false: file cannot be open for writing.
        bool WriteMMCIFData   ( cpstr FName,
                                io::GZ_MODE gzipMode=io::GZM_CHECK );

        /// \brief Writes (next) mmCIF data block into file.
        /// \param RCFile reference to a CFile object opened on a file.
        ///
        /// This function allows for sequential write of mmCIF data blocks
        /// into a file.
        ///
        /// \note Alternatively, files with multiple data blocks can be
        /// created using mmcif::File class.
        ///
        /// Example:
        /*!
      \code
      io::File       f;
      mmcif::Data  cifData;

        // open file first
        f.assign ( "/path/example.cif" );
        if (!f.rewrite())  {
          printf ( " *** cannot open file '%s' for writing.\n",
                   f.FileName() );
          return -1;
        }

        cifData.PutDataName ( "name1" );
        // fill cifData with all data needed
        cifData.WriteMMCIF ( f ); // first data block written

        cifData.FreeMemory  ( 0 );  // reset data block to empty
        cifData.PutDataName ( "name2" );
        // fill cifData with all data needed
        cifData.WriteMMCIF ( f );  // second data block written

        // add as many data blocks as needed

        // now close the file
        f.shut();

      \endcode

        */
        void  WriteMMCIF ( io::RFile f );


        // -------- Retrieving data

        /// \brief Returns the number of categories (structures and loops)
        ///        in data block.
        inline int   GetNumberOfCategories ()  { return nCategories; }

        /// \brief Retrieves pointer to category (a structure or a loop) by
        ///        category number.
        /// \param categoryNo category number to retrieve. Categories are
        ///        numbered from 0 to GetNumberOfCategories()-1.
        /// \return pointer to category, if \b categoryNo is in the right
        ///        range, or \b NULL otherwise.
        ///
        /// \note The category type (structure or loop) is returned by
        /// function mmcif::Category::GetCategoryID().
        /// \note The application should never attempt to deallocate
        /// the category returned. It will be properly disposed of by
        /// mmcif::Data's destructor.
        PCategory GetCategory ( int categoryNo ); // 0..nCategories-1

        /// \brief Retrieves mmCIF structure with given name.
        /// \param CName character string with name of the structure (must
        ///        start with underscore).
        /// \return pointer to structure if structure with given name was
        ///        found, and \b NULL otherwise.
        /// \note The application should never attempt to deallocate
        /// the structure returned. It will be properly disposed of by
        /// mmcif::Data's destructor.
        PStruct GetStructure  ( cpstr CName );

        /// \brief Retrieves mmCIF loop with given name.
        /// \param CName character string with name of the loop (must
        ///        start with underscore).
        /// \return pointer to loop if loop with given name was
        ///        found, and \b NULL otherwise.
        /// \note The application should never attempt to deallocate
        /// the loop returned. It will be properly disposed of by
        /// mmcif::Data's destructor.
        PLoop GetLoop ( cpstr CName );

        /// \brief Finds loop containing all tags from the tag list
        ///        provided.
        /// \param tagList list of tags to be looked for. The list should
        ///        be terminated by empty string \"\". The order of tags
        ///        is not significant.
        /// \return pointer to loop if loop with given tags was found, and
        ///         \b NULL otherwise.
        ///
        /// The function will look for first loop that includes all tags
        /// from the list. The list does not have to include all tags for
        /// that loop in order for function to succeed. This function is
        /// useful for reading \"dirty cifs\" that may contain loops without
        /// a name.
        PLoop FindLoop ( cpstr * tagList );

        /// \brief Retrieves data block name into dynamically-allocated
        ///        string.
        /// \param dname pointer reference to a string that accepts data
        ///        block name. If \b dname is not \b NULL, it is treated
        ///        as a pre-allocated string, which is disposed before
        ///        copying. The application is responsible for deallocating
        ///        \b dname.
        /// \param Remove flag to remove name from the data block.
        void GetDataName ( pstr & dname, bool Remove=false );

        /// \brief Returns data block name.
        inline pstr GetDataName()  { return name; }

        //   CheckData(..) returns positive value if the field is in the
        // file:
        //   CIFRC_Structure  category CName is a structure
        //   CIFRC_Loop       category CName is a loop
        // Negative returns mean:
        //   CIFRC_StructureNoTag  category CName is present,
        //                        it is a structure, but it does not
        //                        have tag TName
        //   CIFRC_LoopNoTag       category CName is present,
        //                        it is a loop, but it does not have
        //                        tag TName
        //   CIFRC_NoCategory      category CName is not present.
        // If TName is set to NULL then only the CName is checked and
        // possible returns are CIFRC_Structure, CIFRC_Loop and
        // CIFRC_NoCategory.
        int  CheckData       ( cpstr CName, cpstr TName );

        int  DeleteCategory  ( cpstr CName );
        int  DeleteStructure ( cpstr CName );
        int  DeleteLoop      ( cpstr CName );

        //   Optimize() optimizes the CIF data in memory allocation. It is
        // a good idea to call it once after extraction of data (GetXXXXXX
        // functions) with Remove flag set on has been completed.
        void Optimize();

        //   GetString(..), GetReal(..) and GetInteger(..) return 0 if the
        // requested field was found and successfully converted. Negative
        // returns mean:
        //    CIFRC_WrongFormat   the field was found but failed to convert
        //                        due to improper numeric format
        //    CIFRC_NoTag         category CName was found, but it does not
        //                        have tag TName
        //    CIFRC_NoCategory    category CName was not found
        //    CIFRC_NotAStructure category CName was found, but it is
        //                        a loop rather than a structure.
        //   GetString(..) will try to dispose Dest unless it is assigned
        // NULL value before the call. The string will be then dynamically
        // allocated and copied.
        //   If Remove is set to true, the field will be removed after
        // extraction.
        int  GetString   ( pstr & Dest, cpstr CName, cpstr TName,
                                        bool Remove=false );
        pstr GetString   ( cpstr CName, cpstr TName, int & RC );
        int  DeleteField ( cpstr CName, cpstr TName );
        int  GetReal     ( realtype & R, cpstr CName,
                           cpstr TName, bool Remove=false );
        int  GetInteger  ( int & I, cpstr CName, cpstr TName,
                                    bool Remove=false );

        //   GetLoopLength(..) returns CIFRC_NotALoop if the category CName
        // is not a loop, CIFRC_NoCategory if the category CName is not
        // found. Non-negative returns give the length of the loop (may be
        // 0 if the loop is empty).
        int  GetLoopLength ( cpstr CName );

        //   GetLoopString(..), GetLoopReal(..) and GetLoopInteger(..) act
        // like GetString(..), GetReal(..) and GetInteger(..) above for
        // nrow-th element of the 'loop_' (indexed like 0..N-1 where N
        // is obtained through GetLoopLength(..)). They will return
        // CIFRC_WrongIndex if nrow is out of range.
        //   If Remove is set to true, the field will be removed after
        // extraction.
        int  GetLoopString   ( pstr & Dest, cpstr CName,
                                            cpstr TName, int nrow,
                                            bool Remove=false );
        pstr GetLoopString   ( cpstr CName, cpstr TName,
                               int nrow, int & RC );
        int  DeleteLoopField ( cpstr CName, cpstr TName,
                               int nrow );
        int  GetLoopReal     ( realtype & R, cpstr CName,
                                             cpstr TName, int nrow,
                                             bool Remove=false );
        int  GetLoopInteger  ( int & I, cpstr CName,
                                        cpstr TName, int nrow,
                                        bool Remove=false );

        //   GetLoopSVector(..), GetLoopRVector(..) and GetLoopIVector(..)
        // read CIF 'loop_' data into allocated vectors of strings, reals
        // and integers, correspondingly. The vectors may be deallocated
        // prior to call and assigned NULL, in which case they will be
        // allocated with offsets of i1, which is also the lower index of
        // the 'loop_' data transferred into it. The upper vector index is
        // given by i2 or by the loop's length whichever is less. If
        // vectors are not assigned NULL prior the call, it is assumed
        // that they are properly (i1-offset, i2-i1+1 length) allocated.
        //   The return codes are same as those of GetLoopString(..),
        // GetLoopReal(..) and GetLoopInteger(..).
        int  GetLoopSVector ( psvector & S, cpstr CName,
                              cpstr TName, int i1=0, int i2=MaxInt4,
                              bool Remove=false );
        int  GetLoopRVector ( rvector  & R, cpstr CName,
                              cpstr TName, int i1=0, int i2=MaxInt4,
                              bool Remove=false );
        int  GetLoopIVector ( ivector  & I, cpstr CName,
                              cpstr TName, int i1=0, int i2=MaxInt4,
                              bool Remove=false );


        // -------- Storing data

        //   Unless the data are to be added to the existing CIF structure,
        // FreeMemory() should be called once before creating a new
        // CIF data set.
        void FreeMemory ( int key );

        void PutDataName ( cpstr dname ); // stores name for 'data_'
                                          // record

        //   PutString(..), PutReal(..) and PutInteger(..) will put the
        // values given into the specified category (CName) under the
        // specified tag (TName). The category, tag and field are created
        // automatically; the field will be replaced silently if identical
        // CName.TName is specified in two calls. Calls of these functions
        // may follow in random order; however CIF file will have all tags
        // grouped by categories and catgories will follow in the order
        // of first appearance in PutString(..), PutReal(..) or
        // PutInteger(..).
        //   Return code - one of CIFRC_Ok or CIFRC_NotAStruct
        int  PutNoData   ( int NoDataType, cpstr CName,
                           cpstr TName );
        int  PutString   ( cpstr S, cpstr CName,
                           cpstr TName, bool Concatenate=false );
        int  PutDate     ( cpstr CName, cpstr TName );
        int  PutReal     ( realtype R, cpstr CName, cpstr TName,
                                       int prec=8 );
        int  PutInteger  ( int I, cpstr CName, cpstr TName );

        //   If loop category CName is not present in the CIF data
        // structure, AddLoop(..) creates an empty one and returns
        // its pointer in Loop. If loop category CName is already in
        // the CIF data structure, its pointer is returned, and any
        // data which might be contained in it, remains untouched.
        //   To stuff the loop with data, first the data tags have to
        // be specified by calling  Loop->AddLoopTag(..). After all
        // tags are given, the data comes as a stream of calls
        // Loop->AddString(..), Loop->AddReal(..) and
        // Loop->AddInteger(..) which should provide data for every
        // tag in sequence in strictly the same order as the tags
        // were given. This essentially reflects reading a CIF loop
        // from a file.
        //   Alternatively, the loop data may be stored with PutLoopXXX()
        // functions given below, although this way may be less
        // efficient (but more flexible).
        //   AddLoop(..) may return
        //     CIFRC_Ok       category was present
        //     CIFRC_Created  category was not present but it has
        //                    been created; the category is empty
        //     CIFRC_NotALoop category was present as a structure, but
        //                    has been replaced for a loop;
        //                    the category is empty.
        int  AddLoop      ( cpstr CName, PLoop   & cifLoop   );
        int  AddStructure ( cpstr CName, PStruct & cifStruct );

        //   PutLoopString(..), PutLoopReal(..) and PutLoopInteger(..) act
        // like PutString(..), PutReal(..) and PutInteger(..) above for
        // nrow-th element of the 'loop_' CName (indexed begining from 0).
        // In consequitive calls, given values of nrow does not have to be
        // ordered; the most efficient way is to start with HIGHEST value
        // for nrow in the loop and move down to 0. The least efficient way
        // is to start with nrow=0 and move up.
        //   These functions allow to form loops in arbitrary way.
        //   The functions may return CIFRC_Ok or CIFRC_NotALoop.
        int  PutLoopNoData  ( int NoDataType, cpstr CName,
                                              cpstr TName, int nrow );
        int  PutLoopString  ( cpstr S,   cpstr CName,
                                              cpstr TName, int nrow );
        int  PutLoopReal    ( realtype R, cpstr CName,
                                          cpstr TName, int nrow,
                                          int  prec=8 );
        int  PutLoopInteger ( int I, cpstr CName, cpstr TName,
                                     int nrow );

        //   PutLoopSVector(..), PutLoopRVector(..) and PutLoopIVector(..)
        // put vectors of values into specified loop fields. Parameters i1
        // and i2 give the range of indices of values which are to be
        // transfered. To transfer an entire vector allocated as [0..N-1]
        // i1 shoudl be set to 0 and i2 - to N-1. Note that the loop is
        // always indexed as starting form 0 on, therefore negative i1 and
        // i2 are not allowed, and specifying i1>0 will leave first i1
        // elements of the CIF loop for the corresponding tag undefined
        // (will be output like '?').
        //   These functions allow to form loops in arbitrary way.
        int  PutLoopSVector ( psvector S, cpstr CName,
                              cpstr TName, int i1, int i2 );
        int  PutLoopRVector ( rvector  R, cpstr CName,
                              cpstr TName, int i1, int i2,
                              int prec=8 );
        int  PutLoopIVector ( ivector  I, cpstr CName,
                              cpstr TName, int i1, int i2 );

        int  RenameCategory ( cpstr CName, cpstr newCName );

        // --------

        void Copy         ( PData Data );
        int  CopyCategory ( PData Data, cpstr CName,
                                        cpstr newCName=NULL );

        void PrintCategories();  // for debuging only

        void write ( io::RFile f );
        void read  ( io::RFile f );

      protected:
        pstr       name;
        int        nCategories;
        PPCategory Category;
        ivector    index;
        int        flags;
        int        Warning;
        int        loopNo;  // used locally for suggesting categories
        int        tagNo;   // used locally for suggesting tags
        psvector   WrongCat;
        psvector   WrongTag;
        int        nWrongFields;

        void  InitData        ();
        void  FreeWrongFields ();
        bool  CheckWrongField ( cpstr C, cpstr T );
        void  Sort            ();

        //   GetCategoryNo searches for index of category cname
        // in Category[]. Return:
        //    >=0 : position of the category found
        //     <0 : the category was not found, it could be inserted before
        //          (-RC-1)th element, where RC is the return value
        int  GetCategoryNo  ( cpstr cname );
        int  AddCategory    ( cpstr cname );
        int  DeleteCategory ( int  CatNo );

        void GetDataItem    ( io::RFile f, pstr S, pstr & L, pstr & p,
                                        int & lcount, int & llen );
        void GetLoop        ( io::RFile f, pstr S, pstr & L, pstr & p,
                                        int & lcount, int & llen );
        int  GetField       ( io::RFile f, pstr S, pstr & L, pstr & p,
                                        int & lcount, int & llen );

    };



    //  ========================  File  =============================

    DefineClass(File);
    DefineStreamFunctions(File);

    class File : public io::Stream  {

      public :
        int      nData;
        ivector  index;
        PPData   data;

        File ();
        File ( cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK );
        File ( io::RPStream Object );
        ~File();

        void  SetPrintWarnings ( bool SPW ) { PrintWarnings = SPW; }
        void  SetStopOnWarning ( bool SOW ) { StopOnWarning = SOW; }

        int   ReadMMCIFFile  ( cpstr FName,
                               io::GZ_MODE gzipMode=io::GZM_CHECK );
        int   WriteMMCIFFile ( cpstr FName,
                               io::GZ_MODE gzipMode=io::GZM_CHECK );

        int   GetNofData()  { return nData; }
        PData GetCIFData     ( int   dataNo );  // 0..nData-1
        PData GetCIFData     ( cpstr DName  );
        int   AddCIFData     ( cpstr DName  );
        int   DeleteCIFData  ( cpstr DName  );
        int   DeleteCIFData  ( int   dataNo );
        int   GetCIFDataNo   ( cpstr DName  );

        void  WriteMMCIF     ( io::RFile f  );

        void  Copy  ( PFile File );

        void  write ( io::RFile f );
        void  read  ( io::RFile f );

      protected:
        int  nAllocData;
        bool PrintWarnings;
        bool StopOnWarning;

        void  InitFile   ();
        void  FreeMemory ();
        void  Sort       ();
        void  ExpandData ( int nDataNew );

    };


    extern MMDB_IMEX pstr GetMMCIFInputBuffer ( int & LineNo );

    //  isCIF will return
    //    -1   if file FName does not exist
    //     0   if file FName is likely a CIF file ( 'data_' is present )
    //     1   if file FName is not a CIF file ( 'data_' is absent )
    extern int isCIF ( cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK );
    extern int isCIF ( io::RFile f );

	MMDB_IMEX pstr GetCIFMessage ( pstr M, int RC );


  }  // namespace mmcif

}  // namespace mmdb


#endif


