#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "pileup.h"
#include "getopt.h"
#define MAX_FN_LEN (512)
#define MAP_QUAL_CUT (30)
#define MIN_ALLELE_BOTH_STRANDS (1) // The minimum number of times we 
// must see an alternate allele on BOTH STRANDS to report a variant position
// This is the default value. Can be changed by an option

int variant( PulP pp, char* allele1, char* allele2,
	     const int het_only, const int mabs );

void help( void ) {
  printf( "pu-variants -p <pileup fn> -l <low cov> -h <high cov> -d \n" );
  printf( "            -m [required observations of each allele on each\n" );
  printf( "                strands; default = %d]\n", MIN_ALLELE_BOTH_STRANDS );
  printf( "            -H [if set, show het sites only]\n" );
  printf( "Prints a bed file of positions with likely variants.\n" );
  printf( "These variants are either homozygous differences from\n" );
  printf( "the reference or heterozygotes within this individual.\n" );
  printf( "To be reported, a variant base must be seen in both\n" );
  printf( "forward and reverse strand at least once.\n" );
  printf( "It is understood that map-quality and base-quality\n" );
  printf( "cut-offs are to be set by the options to samtools mpileup\n" );
  printf( "If the -d flag is set, it runs in debug mode, sending to\n" );
  printf( "stderr the raw mpileup lines that are determined to have\n" );
  printf( "evidence for a variant.\n" );
  exit( 0 );
}
int main( int argc, char* argv[] ) {
  extern char* optarg;
  char pu_fn[MAX_FN_LEN];
  char line[MAX_LINE_LEN+1];
  int ich;
  int pu_file_input = 0;
  int low_cov = 6;
  int high_cov = 20;
  int mabs = MIN_ALLELE_BOTH_STRANDS; // set default
  int debug = 0;
  int het_only = 0;
  int invalid_pul;
  QcutsP qcp;
  PulP pp;
  char allele1, allele2;
  FILE* f;
  char* s1;
  size_t i;

  qcp = dummyQcutsP();

    /* Get options */
  while( (ich=getopt( argc, argv, "p:l:h:m:dH" )) != -1 ) {
    switch(ich) {
    case 'p' :
      strcpy( pu_fn, optarg );
      pu_file_input = 1;
      break;
    case 'l' :
      low_cov = atoi( optarg );
      break;
    case 'h' :
      high_cov = atoi( optarg );
      break;
    case 'm' :
      mabs = atoi( optarg );
      break;
    case 'd' :
      debug = 1;
      break;
    case 'H' :
      het_only = 1;
      break;
    default :
      help();
    }
  }
  if ( argc == 1 ) {
    help();
  }

  /* Initialize pileup lines structures for the input */
  pp = (PulP)malloc(sizeof(Pul));

  /* Decide on the filehandle on the pileup input */
  if ( pu_file_input ) {
    f = fileOpen( pu_fn, "r" );
  }
  else {
    f = stdin;
  }

  s1 = fgets( line, MAX_LINE_LEN, f );
  if ( s1 != NULL ) {
    invalid_pul = line2pul( line, pp );
  }

  /* Now, go through the pu input */
  while( (s1 != NULL) ) {
    if ( !invalid_pul ) {
      if ( (pp->cov <= high_cov) &&
	   (pp->cov >= low_cov) ) {
	if ( variant( pp, &allele1, &allele2, het_only, mabs ) ) {
	  printf( "%s %d %d %c %c\n", 
		  pp->chr, (int)pp->pos, (int)pp->pos,
		  allele1, allele2 );
	  if ( debug ) {
	    fprintf( stderr, "** Pile-up line that has a variant\n%s\n", 
		     line );
	  }
	}
      }
    }
    s1 = fgets( line, MAX_LINE_LEN, f );
    if ( s1 != NULL ) {
      invalid_pul = line2pul( line, pp );
    }
  }


  /* It's a fire hazard - shut it down! */
  fclose( f );
}

/* Takes a PulP pp - structure with the information from
   a pileup line - and returns true if this site has
   (1) at least two alleles with forward and reverse 
   strand instances of both bases 
   or
   (2) a forward and reverse strand instance of an allele
   that is not the reference allele
   If it returns true, it also populates allele1 and allele2
   with the two alleles. allele1 is always the reference
   if that is one of the two alleles
*/
int variant( PulP pp, char* allele1, char* allele2,
	     const int het_only, const int mabs ) {
  int fA = 0;
  int fC = 0;
  int fG = 0;
  int fT = 0;
  int rA = 0;
  int rC = 0;
  int rG = 0;
  int rT = 0;
  int i;
  int num_conf_alleles = 0;
  char conf_alleles[4];

  for( i = 0; i < pp->cov; i++ ) {
    switch (pp->bases[i]) {
    case 'A' :
      if ( pp->strands[i] == 1 ) {
	fA++;
      }
      else {
	rA++;
      }
      break;
    case 'C' :
      if ( pp->strands[i] == 1 ) {
	fC++;
      }
      else {
	rC++;
      }
      break;
    case 'G' :
      if ( pp->strands[i] == 1 ) {
	fG++;
      }
      else {
	rG++;
      }
      break;
    case 'T' :
      if ( pp->strands[i] == 1 ) {
	fT++;
      }
      else {
	rT++;
      }
      break;
    }
  }

  /* Populate the conf_alleles[] and num_conf_alleles */
  if ( (fA >= mabs) && 
       (rA >= mabs) ) {
      conf_alleles[num_conf_alleles] = 'A';
      num_conf_alleles++;
    }
    if ( (fC >= mabs) && 
	 (rC >= mabs) ) {
      conf_alleles[num_conf_alleles] = 'C';
      num_conf_alleles++;
    }
    if ( (fG >= mabs) && 
	 (rG >= mabs) ) {
      conf_alleles[num_conf_alleles] = 'G';
      num_conf_alleles++;
    }
    if ( (fT >= mabs) && 
	 (rT >= mabs) ) {
      conf_alleles[num_conf_alleles] = 'T';
      num_conf_alleles++;
    }
    
    /* Bi-allelic in confident alleles? */
    if ( num_conf_alleles == 2 ) {
      *allele1 = conf_alleles[0];
      *allele2 = conf_alleles[1];
      if ( conf_alleles[0] == pp->ref ) {
	// Ready to return
	return 1;
      }
      if ( conf_alleles[1] == pp->ref ) {
	*allele1 = pp->ref;
	*allele2 = conf_alleles[0];
	return 1;
      }
      /* Neither were the reference, no guaranteed order */
      return 1;
    }

    /* If we only found one confident allele, this is
       a variant site IFF that wasn't the reference
       allele & if we're not in het_only mode */
    if ( het_only ) {
      return 0;
    }
    if ( num_conf_alleles == 1 ) {
      if ( (conf_alleles[0] != pp->ref) &&
	   (pp->ref != 'N') ) {
	*allele1 = pp->ref;
	*allele2 = conf_alleles[0];
	return 1;
      }
    }

    return 0; // tri- or quad-allelic or only reference
}
  
      
