!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc -*- fortran -*-
!
!     SILO include file for Fortan applications.
!
! Programmer:   Robb Matzke
!       Wed Nov 30 11:40:46 PST 1994
!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
!    This file is generated automatically from C header files in the
!    SILO distribution.  Do not make changes to this file manually,
!    as those changes will disappear.
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
!
! This file was built by saying:
!       mv silo.inc silo.inc~
!       cat silo.h silo_f.h | mkinc > silo.inc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
!
   INTERFACE
         FUNCTION dbcreate (NAME, IL, I1, I2, NAME2, IL2, I3, I4)
         !DEC$ ATTRIBUTES ALIAS:'_dbcreate_' :: dbcreate
         INTEGER I1, I2, I3, I4, IL, IL2
         CHARACTER*(*):: NAME, NAME2
         END FUNCTION dbcreate
      END INTERFACE
      
      INTERFACE
         FUNCTION dbmkoptlist (I1, I2)
         !DEC$ ATTRIBUTES ALIAS:'_dbmkoptlist_' :: dbmkoptlist
         INTEGER I1, I2
         END FUNCTION dbmkoptlist
      END INTERFACE
      
      INTERFACE
         FUNCTION dbaddcopt (I1, I2, NAME, IL)
         !DEC$ ATTRIBUTES ALIAS: '_dbaddcopt_' :: dbaddcopt
         INTEGER I1, I2, IL 
         CHARACTER*(*):: NAME
         END FUNCTION dbaddcopt
      END INTERFACE            

      INTERFACE
         FUNCTION dbputqm (I1, NAME1, IL1, NAME2, IL2, NAME3, IL3,    &
     &      NAME4, IL4, X, Y, Z, ND, NN, I2, I3, I4, I5)
         !DEC$ ATTRIBUTES ALIAS:'_dbputqm_' :: dbputqm
         INTEGER I1, I2, I3, I4, I5, NN, IL1, IL2, IL3, IL4
         INTEGER, DIMENSION(*):: ND
         REAL, DIMENSION(*):: X, Y, Z           
         CHARACTER*(*):: NAME1, NAME2, NAME3, NAME4
         END FUNCTION dbputqm
      END INTERFACE
      
      INTERFACE
         FUNCTION dbputqv1 (I1, NAME1, IL1, NAME2, IL2, VAR,    &
     &       ND, NN, MIXVAR, I3, I4, I5, I6, I7)
         !DEC$ ATTRIBUTES ALIAS:'_dbputqv1_' :: dbputqv1
         INTEGER I1, I2, I3, I4, I5, I6, I7, NN, IL1, IL2
         INTEGER, DIMENSION(*):: ND
         REAL, DIMENSION(*):: VAR
         INTEGER, DIMENSION(*):: MIXVAR
         CHARACTER*(*):: NAME1, NAME2
         END FUNCTION dbputqv1
      END INTERFACE      
          
      INTERFACE
         FUNCTION dbfreeoptlist (I)
         !DEC$ ATTRIBUTES ALIAS:'_dbfreeoptlist_':: dbfreeoptlist
         INTEGER I              
         END FUNCTION dbfreeoptlist
      END INTERFACE

      
      INTERFACE
         FUNCTION dbputmat (I1, NAME1, I2, NAME2, I3, NM,                  &
     &       MM, MLIST, MDS, NDS, I4, I5, I6, I7, I8, I9, I10, I11)
         !DEC$ ATTRIBUTES ALIAS:'_dbputmat_' :: dbputmat
         INTEGER I1, I2, I3, I8, I9, I10, I11, NM, NDS
         INTEGER, DIMENSION(*):: MM, MLIST, MDS, I4, I5, I6, I7
         CHARACTER*(*):: NAME1, NAME2
         END FUNCTION dbputmat
      END INTERFACE        
      
      INTERFACE
         FUNCTION dbclose (I)
         !DEC$ ATTRIBUTES ALIAS:'_dbclose_' :: dbclose
         INTEGER I 
         END FUNCTION dbclose
      END INTERFACE
!      
!...Functions (removed: dbcreate, dbclose, dbputqv1, dbputqm, 
!              dbmkoptlist, dbaddcopt, dbfreeoptlist, dbputmat)
!
!      external dbclose, dbfreeoptlist, dbputqv1, dbputqm
!      external dbaddcopt, dbmkoptlist, dbcreate, dbputmat
      external dbaddcaopt, dbaddccmp, dbadddopt, dbaddicmp
      external dbaddiopt, dbaddrcmp, dbaddropt, dbcalcfl 
      external dbfgetca, dbfreeobject 
      external dbget2dstrlen, dbgetca, dbgetcksums, dbgetcompress
      external dbgetcurve, dbgethdfnms, dbgetovrwrt, dbgetqv1, dbinqca
      external dbinqfile, dbinqlen, dbmkdir, dbmkobject
      external dbopen, dbputca, dbputcurve, dbputdefvars, dbputfl
      external dbputmmat, dbputmmesh, dbputmsp, dbputmvar
      external dbputpm, dbputpv1, dbputum, dbputuv1
      external dbputzl, dbputzl2, dbrdvar, dbrdvarslice, dbset2dstrlen
      external dbsetcksums, dbsetcompress, dbsetdir, dbsetdirid
      external dbsethdfnms, dbsetovrwrt, dbshowerrors, dbwrite
      external dbwriteslice, dbwrtfl

!      integer  dbclose, dbfreeoptlist, dbputqv1, dbputqm
!      integer  dbaddcopt, dbmkoptlist, dbcreate, dbputmat
      integer  dbaddcaopt, dbaddccmp, dbadddopt, dbaddicmp
      integer  dbaddiopt, dbaddrcmp, dbaddropt, dbcalcfl 
      integer  dbfgetca, dbfreeobject 
      integer  dbget2dstrlen, dbgetca, dbgetcksums, dbgetcompress
      integer  dbgetcurve, dbgethdfnms, dbgetovrwrt, dbgetqv1, dbinqca
      integer  dbinqfile, dbinqlen, dbmkdir, dbmkobject
      integer  dbopen, dbputca, dbputcurve, dbputdefvars, dbputfl
      integer  dbputmmat, dbputmmesh, dbputmsp, dbputmvar
      integer  dbputpm, dbputpv1, dbputum, dbputuv1
      integer  dbputzl, dbputzl2, dbrdvar, dbrdvarslice, dbset2dstrlen
      integer  dbsetcksums, dbsetcompress, dbsetdir, dbsetdirid
      integer  dbsethdfnms, dbsetovrwrt, dbshowerrors, dbwrite
      integer  dbwriteslice, dbwrtfl

!...Parameters.
      integer  DBCSG_ANGLE_PNLA, DBCSG_ANGLE_PPA, DBCSG_BOX_XYXY
      integer  DBCSG_BOX_XYZXYZ, DBCSG_CIRCLE_PR, DBCSG_COMPLIMENT
      integer  DBCSG_CONE_PNLA, DBCSG_CONE_PPA, DBCSG_CYLINDER_PNLR
      integer  DBCSG_CYLINDER_PPR, DBCSG_DIFF, DBCSG_ELLIPSE_PRR
      integer  DBCSG_ELLIPSOID_PRRR, DBCSG_HEX_6F, DBCSG_INNER
      integer  DBCSG_INTERSECT, DBCSG_LINE_G, DBCSG_LINE_PN
      integer  DBCSG_LINE_PP, DBCSG_LINE_X, DBCSG_LINE_Y, DBCSG_ON
      integer  DBCSG_OUTER, DBCSG_PLANE_G, DBCSG_PLANE_PN
      integer  DBCSG_PLANE_PPP, DBCSG_PLANE_X, DBCSG_PLANE_Y
      integer  DBCSG_PLANE_Z, DBCSG_POLYGON_KP, DBCSG_POLYHEDRON_KF
      integer  DBCSG_PRISM_5F, DBCSG_PYRAMID_5F, DBCSG_QUADRATIC_G
      integer  DBCSG_QUADRIC_G, DBCSG_QUAD_4P, DBCSG_SPHERE_PR
      integer  DBCSG_SWEEP, DBCSG_TET_4F, DBCSG_TRI_3P, DBCSG_UNION
      integer  DBCSG_XFORM, DBOPT_ALIGN, DBOPT_ALLOWMAT0
      integer  DBOPT_ASCII_LABEL, DBOPT_BASEINDEX, DBOPT_BLOCKORIGIN
      integer  DBOPT_BNDNAMES, DBOPT_COORDSYS, DBOPT_CYCLE, DBOPT_DTIME
      integer  DBOPT_EXTENTS, DBOPT_EXTENTS_SIZE, DBOPT_FACETYPE
      integer  DBOPT_GROUPINGNAMES, DBOPT_GROUPINGS
      integer  DBOPT_GROUPINGS_SIZE, DBOPT_GROUPNUM, DBOPT_GROUPORIGIN
      integer  DBOPT_HAS_EXTERNAL_ZONES, DBOPT_HIDE_FROM_GUI
      integer  DBOPT_HI_OFFSET, DBOPT_LABEL, DBOPT_LO_OFFSET
      integer  DBOPT_MAJORORDER, DBOPT_MATCOLORS, DBOPT_MATCOUNTS
      integer  DBOPT_MATLISTS, DBOPT_MATNAME, DBOPT_MATNAMES
      integer  DBOPT_MATNOS, DBOPT_MIXLENS, DBOPT_NGROUPS, DBOPT_NMAT
      integer  DBOPT_NMATNOS, DBOPT_NMATSPEC, DBOPT_NODENUM
      integer  DBOPT_NSPACE, DBOPT_ORIGIN, DBOPT_PHZONELIST
      integer  DBOPT_PLANAR, DBOPT_REFERENCE, DBOPT_REGNAMES, DBOPT_TIME
      integer  DBOPT_TOPO_DIM, DBOPT_UNITS, DBOPT_USESPECMF
      integer  DBOPT_XLABEL, DBOPT_XUNITS, DBOPT_XVARNAME, DBOPT_YLABEL
      integer  DBOPT_YUNITS, DBOPT_YVARNAME, DBOPT_ZLABEL
      integer  DBOPT_ZONECOUNTS, DBOPT_ZONENAMES, DBOPT_ZONENUM
      integer  DBOPT_ZUNITS, DBOPT_ZVARNAME, DB_ABORT, DB_ALL, DB_APPEND
      integer  DB_AREA, DB_ARRAY, DB_BNDCENT, DB_CARTESIAN, DB_CHAR
      integer  DB_CLOBBER, DB_COLLINEAR, DB_COLMAJOR, DB_CRAY
      integer  DB_CSGMESH, DB_CSGVAR, DB_CSGZONELIST, DB_CURVE
      integer  DB_CURVILINEAR, DB_CYLINDRICAL, DB_DEBUG, DB_DEFVARS
      integer  DB_DIR, DB_DOUBLE, DB_EDGELIST, DB_F77NULL, DB_FACECENT
      integer  DB_FACELIST, DB_FLOAT, DB_HDF5, DB_HDF5_CORE
      integer  DB_HDF5_MPIO, DB_HDF5_MPIOP, DB_HDF5_SEC2, DB_HDF5_STDIO
      integer  DB_INT, DB_INTEL, DB_INVALID_OBJECT, DB_LOCAL, DB_LONG
      integer  DB_MATERIAL, DB_MATSPECIES, DB_MULTIBLOCKMESH
      integer  DB_MULTIBLOCKVAR, DB_MULTIMAT, DB_MULTIMATSPECIES
      integer  DB_MULTIMESH, DB_MULTIMESHADJ, DB_MULTIVAR, DB_NETCDF
      integer  DB_NFILES, DB_NFILTERS, DB_NOCLOBBER, DB_NODECENT
      integer  DB_NONCOLLINEAR, DB_NONE, DB_NOTCENT, DB_NOTYPE
      integer  DB_NUMERICAL, DB_OFF, DB_ON, DB_OTHER, DB_PDB
      integer  DB_PHZONELIST, DB_POINTMESH, DB_POINTVAR, DB_QUADMESH
      integer  DB_QUADVAR, DB_QUAD_CURV, DB_QUAD_RECT, DB_READ
      integer  DB_RECTILINEAR, DB_RESUME, DB_ROWMAJOR, DB_RS6000, DB_SGI
      integer  DB_SHORT, DB_SPHERICAL, DB_SUN3, DB_SUN4, DB_SUSPEND
      integer  DB_TAURUS, DB_TOP, DB_UCDMESH, DB_UCDVAR, DB_UNKNOWN
      integer  DB_USERDEF, DB_VARIABLE, DB_VARTYPE_ARRAY
      integer  DB_VARTYPE_LABEL, DB_VARTYPE_MATERIAL, DB_VARTYPE_SCALAR
      integer  DB_VARTYPE_SPECIES, DB_VARTYPE_SYMTENSOR
      integer  DB_VARTYPE_TENSOR, DB_VARTYPE_VECTOR, DB_VOLUME
      integer  DB_ZONECENT, DB_ZONELIST, DB_ZONETYPE_BEAM
      integer  DB_ZONETYPE_HEX, DB_ZONETYPE_POLYGON
      integer  DB_ZONETYPE_POLYHEDRON, DB_ZONETYPE_PRISM
      integer  DB_ZONETYPE_PYRAMID, DB_ZONETYPE_QUAD, DB_ZONETYPE_TET
      integer  DB_ZONETYPE_TRIANGLE, SILO_VERS_MAJ, SILO_VERS_MIN
      integer  SILO_VERS_PAT, SILO_VERS_PRE
      character*12  DB_F77NULLSTRING

      parameter (SILO_VERS_MAJ             = 4              )
      parameter (SILO_VERS_MIN             = 6              )
      parameter (SILO_VERS_PAT             = 0              )
      parameter (SILO_VERS_PRE             = 3              )
      parameter (DB_NETCDF                 = 0              )
      parameter (DB_PDB                    = 2              )
      parameter (DB_TAURUS                 = 3              )
      parameter (DB_UNKNOWN                = 5              )
      parameter (DB_DEBUG                  = 6              )
      parameter (DB_HDF5                   = 7              )
      parameter (DB_HDF5_SEC2              = 256            )
      parameter (DB_HDF5_STDIO             = 512            )
      parameter (DB_HDF5_CORE              = 768            )
      parameter (DB_HDF5_MPIO              = 1024           )
      parameter (DB_HDF5_MPIOP             = 1280           )
      parameter (DB_NFILES                 = 256            )
      parameter (DB_NFILTERS               = 32             )
      parameter (DB_INVALID_OBJECT         = -1             )
      parameter (DB_QUADMESH               = 500            )
      parameter (DB_QUADVAR                = 501            )
      parameter (DB_UCDMESH                = 510            )
      parameter (DB_UCDVAR                 = 511            )
      parameter (DB_MULTIMESH              = 520            )
      parameter (DB_MULTIVAR               = 521            )
      parameter (DB_MULTIMAT               = 522            )
      parameter (DB_MULTIMATSPECIES        = 523            )
      parameter (DB_MULTIBLOCKMESH         = DB_MULTIMESH   )
      parameter (DB_MULTIBLOCKVAR          = DB_MULTIVAR    )
      parameter (DB_MULTIMESHADJ           = 524            )
      parameter (DB_MATERIAL               = 530            )
      parameter (DB_MATSPECIES             = 531            )
      parameter (DB_FACELIST               = 550            )
      parameter (DB_ZONELIST               = 551            )
      parameter (DB_EDGELIST               = 552            )
      parameter (DB_PHZONELIST             = 553            )
      parameter (DB_CSGZONELIST            = 554            )
      parameter (DB_CSGMESH                = 555            )
      parameter (DB_CSGVAR                 = 556            )
      parameter (DB_CURVE                  = 560            )
      parameter (DB_DEFVARS                = 565            )
      parameter (DB_POINTMESH              = 570            )
      parameter (DB_POINTVAR               = 571            )
      parameter (DB_ARRAY                  = 580            )
      parameter (DB_DIR                    = 600            )
      parameter (DB_VARIABLE               = 610            )
      parameter (DB_USERDEF                = 700            )
      parameter (DB_INT                    = 16             )
      parameter (DB_SHORT                  = 17             )
      parameter (DB_LONG                   = 18             )
      parameter (DB_FLOAT                  = 19             )
      parameter (DB_DOUBLE                 = 20             )
      parameter (DB_CHAR                   = 21             )
      parameter (DB_NOTYPE                 = 25             )
      parameter (DB_CLOBBER                = 0              )
      parameter (DB_NOCLOBBER              = 1              )
      parameter (DB_READ                   = 1              )
      parameter (DB_APPEND                 = 2              )
      parameter (DB_LOCAL                  = 0              )
      parameter (DB_SUN3                   = 10             )
      parameter (DB_SUN4                   = 11             )
      parameter (DB_SGI                    = 12             )
      parameter (DB_RS6000                 = 13             )
      parameter (DB_CRAY                   = 14             )
      parameter (DB_INTEL                  = 15             )
      parameter (DBOPT_ALIGN               = 260            )
      parameter (DBOPT_COORDSYS            = 262            )
      parameter (DBOPT_CYCLE               = 263            )
      parameter (DBOPT_FACETYPE            = 264            )
      parameter (DBOPT_HI_OFFSET           = 265            )
      parameter (DBOPT_LO_OFFSET           = 266            )
      parameter (DBOPT_LABEL               = 267            )
      parameter (DBOPT_XLABEL              = 268            )
      parameter (DBOPT_YLABEL              = 269            )
      parameter (DBOPT_ZLABEL              = 270            )
      parameter (DBOPT_MAJORORDER          = 271            )
      parameter (DBOPT_NSPACE              = 272            )
      parameter (DBOPT_ORIGIN              = 273            )
      parameter (DBOPT_PLANAR              = 274            )
      parameter (DBOPT_TIME                = 275            )
      parameter (DBOPT_UNITS               = 276            )
      parameter (DBOPT_XUNITS              = 277            )
      parameter (DBOPT_YUNITS              = 278            )
      parameter (DBOPT_ZUNITS              = 279            )
      parameter (DBOPT_DTIME               = 280            )
      parameter (DBOPT_USESPECMF           = 281            )
      parameter (DBOPT_XVARNAME            = 282            )
      parameter (DBOPT_YVARNAME            = 283            )
      parameter (DBOPT_ZVARNAME            = 284            )
      parameter (DBOPT_ASCII_LABEL         = 285            )
      parameter (DBOPT_MATNOS              = 286            )
      parameter (DBOPT_NMATNOS             = 287            )
      parameter (DBOPT_MATNAME             = 288            )
      parameter (DBOPT_NMAT                = 289            )
      parameter (DBOPT_NMATSPEC            = 290            )
      parameter (DBOPT_BASEINDEX           = 291            )
      parameter (DBOPT_ZONENUM             = 292            )
      parameter (DBOPT_NODENUM             = 293            )
      parameter (DBOPT_BLOCKORIGIN         = 294            )
      parameter (DBOPT_GROUPNUM            = 295            )
      parameter (DBOPT_GROUPORIGIN         = 296            )
      parameter (DBOPT_NGROUPS             = 297            )
      parameter (DBOPT_MATNAMES            = 298            )
      parameter (DBOPT_EXTENTS_SIZE        = 299            )
      parameter (DBOPT_EXTENTS             = 300            )
      parameter (DBOPT_MATCOUNTS           = 301            )
      parameter (DBOPT_MATLISTS            = 302            )
      parameter (DBOPT_MIXLENS             = 303            )
      parameter (DBOPT_ZONECOUNTS          = 304            )
      parameter (DBOPT_HAS_EXTERNAL_ZONES  = 305            )
      parameter (DBOPT_PHZONELIST          = 306            )
      parameter (DBOPT_MATCOLORS           = 307            )
      parameter (DBOPT_BNDNAMES            = 308            )
      parameter (DBOPT_REGNAMES            = 309            )
      parameter (DBOPT_ZONENAMES           = 310            )
      parameter (DBOPT_HIDE_FROM_GUI       = 311            )
      parameter (DBOPT_TOPO_DIM            = 312            )
      parameter (DBOPT_REFERENCE           = 313            )
      parameter (DBOPT_GROUPINGS_SIZE      = 314            )
      parameter (DBOPT_GROUPINGS           = 315            )
      parameter (DBOPT_GROUPINGNAMES       = 316            )
      parameter (DBOPT_ALLOWMAT0           = 317            )
      parameter (DB_TOP                    = 0              )
      parameter (DB_NONE                   = 1              )
      parameter (DB_ALL                    = 2              )
      parameter (DB_ABORT                  = 3              )
      parameter (DB_SUSPEND                = 4              )
      parameter (DB_RESUME                 = 5              )
      parameter (DB_ROWMAJOR               = 0              )
      parameter (DB_COLMAJOR               = 1              )
      parameter (DB_COLLINEAR              = 130            )
      parameter (DB_NONCOLLINEAR           = 131            )
      parameter (DB_QUAD_RECT              = DB_COLLINEAR   )
      parameter (DB_QUAD_CURV              = DB_NONCOLLINEAR)
      parameter (DB_NOTCENT                = 0              )
      parameter (DB_NODECENT               = 110            )
      parameter (DB_ZONECENT               = 111            )
      parameter (DB_FACECENT               = 112            )
      parameter (DB_BNDCENT                = 113            )
      parameter (DB_CARTESIAN              = 120            )
      parameter (DB_CYLINDRICAL            = 121            )
      parameter (DB_SPHERICAL              = 122            )
      parameter (DB_NUMERICAL              = 123            )
      parameter (DB_OTHER                  = 124            )
      parameter (DB_RECTILINEAR            = 100            )
      parameter (DB_CURVILINEAR            = 101            )
      parameter (DB_AREA                   = 140            )
      parameter (DB_VOLUME                 = 141            )
      parameter (DB_ON                     = 1000           )
      parameter (DB_OFF                    = -1000          )
      parameter (DB_VARTYPE_SCALAR         = 200            )
      parameter (DB_VARTYPE_VECTOR         = 201            )
      parameter (DB_VARTYPE_TENSOR         = 202            )
      parameter (DB_VARTYPE_SYMTENSOR      = 203            )
      parameter (DB_VARTYPE_ARRAY          = 204            )
      parameter (DB_VARTYPE_MATERIAL       = 205            )
      parameter (DB_VARTYPE_SPECIES        = 206            )
      parameter (DB_VARTYPE_LABEL          = 207            )
      parameter (DBCSG_QUADRIC_G           = 16777216       )
      parameter (DBCSG_SPHERE_PR           = 33619968       )
      parameter (DBCSG_ELLIPSOID_PRRR      = 33685504       )
      parameter (DBCSG_PLANE_G             = 50331648       )
      parameter (DBCSG_PLANE_X             = 50397184       )
      parameter (DBCSG_PLANE_Y             = 50462720       )
      parameter (DBCSG_PLANE_Z             = 50528256       )
      parameter (DBCSG_PLANE_PN            = 50593792       )
      parameter (DBCSG_PLANE_PPP           = 50659328       )
      parameter (DBCSG_CYLINDER_PNLR       = 67108864       )
      parameter (DBCSG_CYLINDER_PPR        = 67174400       )
      parameter (DBCSG_BOX_XYZXYZ          = 83886080       )
      parameter (DBCSG_CONE_PNLA           = 100663296      )
      parameter (DBCSG_CONE_PPA            = 100728832      )
      parameter (DBCSG_POLYHEDRON_KF       = 117440512      )
      parameter (DBCSG_HEX_6F              = 117506048      )
      parameter (DBCSG_TET_4F              = 117571584      )
      parameter (DBCSG_PYRAMID_5F          = 117637120      )
      parameter (DBCSG_PRISM_5F            = 117702656      )
      parameter (DBCSG_QUADRATIC_G         = 134217728      )
      parameter (DBCSG_CIRCLE_PR           = 150994944      )
      parameter (DBCSG_ELLIPSE_PRR         = 151060480      )
      parameter (DBCSG_LINE_G              = 167772160      )
      parameter (DBCSG_LINE_X              = 167837696      )
      parameter (DBCSG_LINE_Y              = 167903232      )
      parameter (DBCSG_LINE_PN             = 167968768      )
      parameter (DBCSG_LINE_PP             = 168034304      )
      parameter (DBCSG_BOX_XYXY            = 184549376      )
      parameter (DBCSG_ANGLE_PNLA          = 201326592      )
      parameter (DBCSG_ANGLE_PPA           = 201392128      )
      parameter (DBCSG_POLYGON_KP          = 218103808      )
      parameter (DBCSG_TRI_3P              = 218169344      )
      parameter (DBCSG_QUAD_4P             = 218234880      )
      parameter (DBCSG_INNER               = 2130706432     )
      parameter (DBCSG_OUTER               = 2130771968     )
      parameter (DBCSG_ON                  = 2130837504     )
      parameter (DBCSG_UNION               = 2130903040     )
      parameter (DBCSG_INTERSECT           = 2130968576     )
      parameter (DBCSG_DIFF                = 2131034112     )
      parameter (DBCSG_COMPLIMENT          = 2131099648     )
      parameter (DBCSG_XFORM               = 2131165184     )
      parameter (DBCSG_SWEEP               = 2131230720     )
      parameter (DB_F77NULL                = (-99)          )
      parameter (DB_ZONETYPE_BEAM          = 10             )
      parameter (DB_ZONETYPE_POLYGON       = 20             )
      parameter (DB_ZONETYPE_TRIANGLE      = 23             )
      parameter (DB_ZONETYPE_QUAD          = 24             )
      parameter (DB_ZONETYPE_POLYHEDRON    = 30             )
      parameter (DB_ZONETYPE_TET           = 34             )
      parameter (DB_ZONETYPE_PYRAMID       = 35             )
      parameter (DB_ZONETYPE_PRISM         = 36             )
      parameter (DB_ZONETYPE_HEX           = 38             )
      parameter (DB_F77NULLSTRING          = "NULLSTRING"   )


! End.