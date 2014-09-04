#ifndef CLASP_H
#define CLASP_H

/**
 * clasp.h
 * fast fragment chaining
 * using sop gap costs
 * 
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Nov  2 10:06:11 CET 2009
 */

/*
 * SVN
 * Revision of last commit: $Rev: 116 $
 * Author: $Author: steve $
 * Date: $Date: 2010-06-30 13:51:27 +0200 (Wed, 30 Jun 2010) $
 * Id: $Id: clasp.h 116 2010-06-30 11:51:27Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/clasp.h $
 */

#include <stdio.h>
#include <stdlib.h>

#define SOP 	((unsigned char) (0 << 0))
#define LIN	((unsigned char) (1 << 0))
#define VERSION "1.1"

/* Typedef */
typedef struct {
  char *infilename;
  char *outfilename;
  FILE *dev;
  Container *fragments;
  Container *lines;
  Container *subject;
  unsigned char chainmode;
  double lambda;
  double epsilon;
  double minscore;
  int maxgap;
  Uint minfrag;
  Uint colnum;
  Uint* colorder;
  Uint* idcol;
  int idcolnum;
  BOOL outputc;
  BOOL outputf;
  BOOL outputm;
  BOOL outputorig;
} claspinfo_t;


inline static void
bl_claspinfoInit(claspinfo_t *info){
  info->infilename = NULL;
  info->outfilename = NULL;
  info->dev = stdout;
  info->fragments = NULL;
  info->lines = NULL;
  info->subject = NULL;
  info->chainmode = SOP;
  info->lambda = 1;           // chaining parameter
  info->epsilon = 0;          // chaining parameter
  info->maxgap = -1;
  info->minscore = 0;
  info->minfrag = 0;
  info->colorder = NULL;
  info->idcol = NULL;
  info->idcolnum = 0;
  info->outputc = 1;
  info->outputf = 0;
  info->outputm = 0;  
  info->outputorig = 0;
}

inline static void
bl_claspinfoDestruct(claspinfo_t *info){
  Uint i;
  if (info->colorder){
    free(info->colorder);
  }
  if (info->idcol){
    free(info->idcol);
  }
  if (info->outfilename){
    fclose(info->dev);
  }
  if (info->fragments){
    for (i = 0; i < bl_containerSize(info->fragments); i++){
      slmatch_t *sl = (slmatch_t *) bl_containerGet(info->fragments, i);
      if (sl->chain != NULL){
	DBG("still at least one chain not freed before end: %d", i);
	exit(-1);
        bl_slchainDestruct(sl->chain);
        free(sl->chain);
      }
    }
    bl_containerDestruct(info->fragments, bl_slmatchDestruct);
    free(info->fragments);
  }
  if (info->lines){
    for (i = 0; i < bl_containerSize(info->lines); i++){
      free(*(char **) bl_containerGet(info->lines, i));
    }
    bl_containerDestruct(info->lines, NULL);
    free(info->lines);
  }
  if (info->subject){
    for (i = 0; i < bl_containerSize(info->subject); i++){
      free(*(char **) bl_containerGet(info->subject, i));
    }
    bl_containerDestruct(info->subject, NULL);
    free(info->subject);
  }
}

/*
void bl_slWriteHeader();
void bl_slmatchInitFromFile(Container *fragments, Container *lines,
			    Container *subject, char *filename,
			    char *delim, Uint *colorder,
			    Uint *idcol, int idcolnum);
*/
#endif /* CLASP_H */
