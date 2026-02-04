#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <float.h> 
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpzi.h>
#include <flint/fmpq.h>
#include <flint/fmpz_factor.h>
#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"
//clear && gcc DescentSieve.c -o m.o -lm -lgmp -lmpfr -lflint && ./m.o 
typedef struct descent_sieve_struct *DescentSieve;
typedef struct descent_factor_struct *DescentFactor;
struct descent_sieve_struct
{
	int factorBaseBitLength;
	fmpz_t t;
	fmpz_t tPrime;
	fmpz_t t_tPrime;
	fmpz_t t_tprime_minp;
	fmpz_t yF2;
	fmpz_t prime;
	fmpz_t rootPrime;
	fmpz_t bestNum;
	fmpz_t bestDen;
	fmpz_t exp2Value;
	double descentCeilBits;
	double descentFloorBits;
	int exp2;
	uint64_t *primes;
	uint64_t *coefficientF1;
	uint64_t *coefficientF2;
	double bestMSE;
	uint64_t bestX;
	uint64_t bestY;
	int bestF1Sign;
	int bestF2Sign;
};

struct descent_factor_struct
{
	fmpz_t p;
	int exponent;
};

DescentSieve DescentSieve_Create(int factorBaseBitLength, fmpz_t prime)
{
	DescentSieve sieve= malloc(sizeof(struct descent_sieve_struct));
	sieve->factorBaseBitLength = factorBaseBitLength;
	fmpz_init(sieve->t);fmpz_init(sieve->tPrime);
	fmpz_init(sieve->t_tprime_minp);fmpz_init(sieve->yF2);
	fmpz_init(sieve->rootPrime);fmpz_init(sieve->prime);
	fmpz_init(sieve->bestNum);fmpz_init(sieve->bestDen);
	fmpz_init(sieve->t_tPrime);
	fmpz_init(sieve->exp2Value);
	
	sieve->primes = NULL;
	sieve->coefficientF1 = NULL;
	sieve->coefficientF2 = NULL;
	fmpz_set(sieve->prime, prime);
	fmpz_sqrt(sieve->rootPrime, prime);
	sieve->descentCeilBits = fmpz_dlog(sieve->rootPrime) / log(2.0);
	sieve->descentFloorBits = sieve->descentCeilBits - 1;
	sieve->exp2 = 0;
	sieve->bestF1Sign = 0;
	sieve->bestF2Sign = 0;
	//Generate primes with libgmp
	mpz_t currentPrime;
	mpz_init(currentPrime);
	mpz_set_ui(currentPrime, 2);
	size_t factorBaseSize = mpz_sizeinbase(currentPrime,2);
	while(factorBaseSize < factorBaseBitLength)
	{
		uint64_t primeno = mpz_get_ui(currentPrime);
		arrput(sieve->primes, primeno);
		mpz_nextprime(currentPrime, currentPrime);
		factorBaseSize = mpz_sizeinbase(currentPrime,2);
	}
	//Init sieve arrays
	sieve->coefficientF1 = calloc(arrlen(sieve->primes), sizeof(uint64_t));
	sieve->coefficientF2 = calloc(arrlen(sieve->primes), sizeof(uint64_t));
	sieve->bestMSE = DBL_MAX;
	sieve->bestX = 0;
	sieve->bestY = 0;
	mpz_clear(currentPrime);
	return sieve;
}

DescentFactor DescentFactor_CreateDescentFactor()
{
	DescentFactor descentFactor = malloc(sizeof(struct descent_factor_struct));
	fmpz_init(descentFactor->p);
	descentFactor->exponent = 0;
	return descentFactor;
}

int DescentFactor_cmp_fmpz(const void *a, const void *b)
{
	const DescentFactor da = *(const DescentFactor *)a;
	const DescentFactor db = *(const DescentFactor *)b;

	return fmpz_cmp(da->p, db->p);
}

void DescentFactor_Destroy(DescentFactor descentFactor)
{
	if(descentFactor)
	{
		fmpz_clear(descentFactor->p);
		free(descentFactor);
	}
}

void DescentSieve_Destroy(DescentSieve sieve)
{
	if(sieve)
	{
		fmpz_clear(sieve->exp2Value);
		fmpz_clear(sieve->t_tPrime);
		fmpz_clear(sieve->bestNum);fmpz_clear(sieve->bestDen);
		fmpz_clear(sieve->rootPrime);fmpz_clear(sieve->prime);
		fmpz_clear(sieve->t);fmpz_clear(sieve->tPrime);
		fmpz_clear(sieve->t_tprime_minp);fmpz_clear(sieve->yF2);
		arrfree(sieve->primes);
		free(sieve->coefficientF1);
		free(sieve->coefficientF2);
		free(sieve);
	}
}

	
void DescentSieve_PrintSieve(DescentSieve descentSieve)
{
	//printf("Prime: ");fmpz_print(descentSieve->prime);printf("\n");
	printf("t: ");fmpz_print(descentSieve->t);printf("\n");
	printf("tPrime: ");fmpz_print(descentSieve->tPrime);printf("\n");
	printf("t_tprime_minp: ");fmpz_print(descentSieve->t_tprime_minp);printf("\n");		
}

void DescentSieve_PrintFactors(fmpz_factor_t factors)
{
	for(slong i = 0; i < factors->num; i++)
	{
		int size = fmpz_sizeinbase(&factors->p[i], 2);
		printf("(|%d| ", size);
		fmpz_print(&factors->p[i]); printf(" ^ ");
		fmpz_print(&factors->exp[i]); printf("), ");
	}
	printf("\n");
}

void DescentSieve_UpdateSieveConstants(DescentSieve descentSieve, fmpz_t target)
{
	fmpz_t temp0,temp1,coefficientF1,coefficientF2;
	fmpz_init(temp0);
	fmpz_init(temp1);
	fmpz_init(coefficientF1);
	fmpz_init(coefficientF2);
	
	//Assumes t is always positive and smaller than square root of prime
	double tLog_Base2 = fmpz_dlog(target) / log(2.0);
	assert(tLog_Base2 <= descentSieve->descentCeilBits);
	
	//Find 2 multiplier
	descentSieve->exp2 = (int) ceil(descentSieve->descentFloorBits - tLog_Base2);
	fmpz_one_2exp(descentSieve->exp2Value, descentSieve->exp2);
	//Find multiplied t
	fmpz_set(descentSieve->t, target);for(int i = 0; i < descentSieve->exp2; i++){fmpz_mul_ui(descentSieve->t, descentSieve->t, 2);}
	//Find t prime
	fmpz_cdiv_q(descentSieve->tPrime, descentSieve->prime, descentSieve->t);
	//Find t_tprime
	fmpz_mul(descentSieve->t_tPrime,descentSieve->t,descentSieve->tPrime);
	fmpz_mod(descentSieve->t_tPrime,descentSieve->t_tPrime,descentSieve->prime);
	//Find t_tprime_minp	
	fmpz_mul(descentSieve->t_tprime_minp, descentSieve->t, descentSieve->tPrime);
	fmpz_sub(descentSieve->t_tprime_minp, descentSieve->t_tprime_minp, descentSieve->prime);
	
	//Set sieve coefficients
	//F1 sieve coefficient is -tPrime
	fmpz_mul_si(coefficientF1, descentSieve->tPrime, -1);
	//F2 sieve coefficient is -t_tprime_minp * tInverse
	fmpz_mul_si(coefficientF2, descentSieve->t_tprime_minp, -1);	
	for(size_t i = 0; i < arrlen(descentSieve->primes); i++)
	{
		fmpz_set_ui(temp1, descentSieve->primes[i]);
		if(fmpz_invmod(temp1, descentSieve->t, temp1))
		{
			fmpz_mul(temp1, temp1, coefficientF2);	
			descentSieve->coefficientF1[i] = fmpz_mod_ui(temp0, coefficientF1, descentSieve->primes[i]);			
			descentSieve->coefficientF2[i] = fmpz_mod_ui(temp1, temp1, descentSieve->primes[i]);						
		}
		else
		{
			//Set to -1
			descentSieve->coefficientF1[i] = UINT64_MAX;
			descentSieve->coefficientF2[i] = UINT64_MAX;
		}
		//printf("%ld: (%lu), %lu,%lu\n", i,descentSieve->primes[i],descentSieve->coefficientF1[i],descentSieve->coefficientF2[i]);	
	}	
	//Reset best
	descentSieve->bestMSE = DBL_MAX;
	descentSieve->bestX = 0;
	descentSieve->bestY = 0;
	fmpz_clear(temp0);
	fmpz_clear(temp0);
	fmpz_clear(temp1);
	fmpz_clear(coefficientF1);
	fmpz_clear(coefficientF2);
}

double DescentSieve_FindCumulativeMSE(int64_t diffF1, int64_t diffF2)
{
	double sum = (double)(diffF1* diffF1) + (double)(diffF2 * diffF2);
	return sum / 2;
}


void DescentSieve_SievePair(DescentSieve descentSieve, uint64_t x, uint64_t y, uint64_t *targetF1, uint64_t *targetF2)
{
	//Reset sieve arrays
	memset(targetF1, 0, x * sizeof(uint64_t));		
	memset(targetF2, 0, x * sizeof(uint64_t));	
	for(size_t i = 0; i < arrlen(descentSieve->primes); i++)
	{
		//Only work with values where t inverse exists
		if(descentSieve->coefficientF1[i] != UINT64_MAX)
		{
			//Update sieve value by this much
			int sieveValue = 1000 * log((double)descentSieve->primes[i]);
			//Find update index mod p
			uint64_t updateIndexF1 = (descentSieve->coefficientF1[i] * y) % descentSieve->primes[i];
			uint64_t updateIndexF2 = (descentSieve->coefficientF2[i] * y) % descentSieve->primes[i];
			//printf("%ld: (%lu), (F1: %lu, %lu), (F2: %lu, %lu),\n", i,descentSieve->primes[i],descentSieve->coefficientF1[i],updateIndexF1, descentSieve->coefficientF2[i], updateIndexF2);	
			//Update sieve arrays
			for(uint64_t sieveIndex = updateIndexF1; sieveIndex < x; sieveIndex += descentSieve->primes[i])
			{
				targetF1[sieveIndex] += sieveValue;
			}
			for(uint64_t sieveIndex = updateIndexF2; sieveIndex < x; sieveIndex += descentSieve->primes[i])
			{
				targetF2[sieveIndex] += sieveValue;
			}
		}
	}
	
	//Find out which is closest to target log
	fmpz_t intF1,intF2,temp0;
	fmpz_init(intF1);fmpz_init(intF2);fmpz_init(temp0);
	
	for(int targetX = 0; targetX < x; targetX++)
	{
		//intF1 = targetX + tPrime * y
		fmpz_mul_ui(intF1, descentSieve->tPrime, y);
		fmpz_add_ui(intF1, intF1, targetX);	
		fmpz_mod(intF1,intF1,descentSieve->prime);
		uint64_t intF1Log = 1000 * (fmpz_dlog(intF1));
		uint64_t sieveF1Log = targetF1[targetX];
		int64_t differenceF1 = (int64_t)intF1Log - (int64_t)sieveF1Log;
		//intF2 = t * targetX + (tt' - p) * y
		fmpz_mul_ui(intF2, descentSieve->t_tprime_minp, y);
		fmpz_mul_ui(temp0, descentSieve->t, targetX);
		fmpz_add(intF2, intF2, temp0);
		fmpz_mod(intF2,intF2,descentSieve->prime);
		uint64_t intF2Log = 1000 * (fmpz_dlog(intF2));
		uint64_t sieveF2Log = targetF2[targetX];
		int64_t differenceF2 = (int64_t)intF2Log - (int64_t)sieveF2Log;
		
		double mse = DescentSieve_FindCumulativeMSE(differenceF1, differenceF2);
		if(mse < descentSieve->bestMSE)
		{
			descentSieve->bestMSE = mse;
			descentSieve->bestX   = targetX;
			descentSieve->bestY   = y;
			
			fmpz_set(descentSieve->bestNum, intF2);
			fmpz_set(descentSieve->bestDen, intF1);
			//printf("MSE: %.3f\nintF1: (%lu, %lu) ",mse,intF1Log, sieveF1Log);fmpz_print(intF1);printf("\n");	
			printf("best MSE: %.3f, X: %lu, Y: %lu\n",descentSieve->bestMSE,descentSieve->bestX,descentSieve->bestY);
		}
		//printf("MSE: %.3f\nintF1: (%lu, %lu) ",mse,intF1Log, sieveF1Log);fmpz_print(intF1);printf("\n");	
		//printf("intF2: (%lu, %lu) ",intF2Log, sieveF2Log);fmpz_print(intF2);printf("\n\n");		
	}
	fmpz_clear(intF1);
	fmpz_clear(intF2);
	fmpz_clear(temp0);
}


DescentFactor *DescentSieve_Factors(fmpz_factor_t factors, DescentSieve descentSieve,uint64_t sieveX,uint64_t sieveY,uint64_t maxYoffset)
{	
	DescentFactor *descentFactors = NULL;
	fmpz_t temp0;
	fmpz_init(temp0);
	fmpz_factor_t factorNum;fmpz_factor_init(factorNum);
	fmpz_factor_t factorDen;fmpz_factor_init(factorDen);
	for(slong i = 0; i < factors->num; i++)
	{
		int size = fmpz_sizeinbase(&factors->p[i], 2);
		bool reductionOccurred = false;
		//printf("(|%d| ", size);fmpz_print(&factors->p[i]);printf("\n");				
		if(size > descentSieve->factorBaseBitLength)
		{
			DescentSieve_UpdateSieveConstants(descentSieve, &factors->p[i]);	
			
			for(int yOffset = 0; yOffset < maxYoffset; yOffset++)			
			{
				//We can try diff y's
				uint64_t *targetF1 = calloc(sieveX, sizeof(uint64_t));
				uint64_t *targetF2 = calloc(sieveX, sizeof(uint64_t));
				DescentSieve_SievePair(descentSieve, sieveX, sieveY+yOffset,targetF1,targetF2);
				//printf("best MSE: %.3f, X: %lu, Y: %lu\n",descentSieve->bestMSE,descentSieve->bestX,descentSieve->bestY);
				//ContradictionTest(descentSieve);
				free(targetF1);
				free(targetF2);
			}
			//Remember to multiply denominator by 2 multipler
			fmpz_mul(temp0, descentSieve->bestDen, descentSieve->exp2Value);
			fmpz_invmod(temp0,temp0, descentSieve->prime);
			fmpz_mul(temp0, temp0, descentSieve->bestNum);
			fmpz_mod(temp0, temp0, descentSieve->prime);		
			assert(fmpz_cmp(temp0, &factors->p[i]) == 0);
			
			fmpz_factor(factorNum, descentSieve->bestNum);
			fmpz_factor(factorDen, descentSieve->bestDen);	
			size_t largestFactor = 0;
			//Quick check to test actual reduction
			for(slong j = 0; j < factorNum->num; j++){size_t currentSize = fmpz_sizeinbase(&factorNum->p[j],2);if(currentSize > largestFactor){largestFactor = currentSize;}}
			for(slong j = 0; j < factorDen->num; j++){size_t currentSize = fmpz_sizeinbase(&factorDen->p[j],2);if(currentSize > largestFactor){largestFactor = currentSize;}}
			
			if(largestFactor > size)
			{
				//printf("ggFyck\n");
				reductionOccurred = false;
			}
			else
			{
				reductionOccurred = true;
				//Create DescentFactors
				DescentFactor twoDenominator = DescentFactor_CreateDescentFactor();
				fmpz_set_ui(twoDenominator->p, 2);
				twoDenominator->exponent = (int)descentSieve->exp2;	
				twoDenominator->exponent *= -1;
				arrput(descentFactors, twoDenominator);
				for(slong j = 0; j < factorNum->num; j++)
				{
					DescentFactor descentFactor = DescentFactor_CreateDescentFactor();
					fmpz_set(descentFactor->p, &factorNum->p[j]);
					descentFactor->exponent = factorNum->exp[j];	
					arrput(descentFactors, descentFactor);
				}
				
				for(slong j = 0; j < factorDen->num; j++)
				{
					DescentFactor descentFactor = DescentFactor_CreateDescentFactor();
					fmpz_set(descentFactor->p, &factorDen->p[j]);
					descentFactor->exponent = (int)factorDen->exp[j];	
					descentFactor->exponent *= -1;
					arrput(descentFactors, descentFactor);
				}
				
				
				//printf("(|%d| ", size);fmpz_print(&factors->p[i]); printf(" ^ ");DescentSieve_PrintSieve(descentSieve);
				//printf("BestNum: ");fmpz_print(descentSieve->bestNum);printf("\n");
				//printf("BestDen: ");fmpz_print(descentSieve->bestDen);printf("\n");
				//printf("Res: ");fmpz_print(temp0);printf("\n");
			}
			//printf("best MSE: %.3f, X: %lu, Y: %lu\n",descentSieve->bestMSE,descentSieve->bestX,descentSieve->bestY);
			//ContradictionTest(descentSieve);
		}
		
		if(reductionOccurred == false)
		{
			DescentFactor descentFactor = DescentFactor_CreateDescentFactor();
			fmpz_set(descentFactor->p, &factors->p[i]);
			descentFactor->exponent = factors->exp[i];	
			arrput(descentFactors, descentFactor);
		}	
	}
	//printf("\n");
	
	fmpz_clear(temp0);
	fmpz_factor_clear(factorNum);fmpz_factor_clear(factorDen);
	qsort(descentFactors,arrlen(descentFactors),sizeof(DescentFactor),DescentFactor_cmp_fmpz);
	//for(size_t i = 0; i < arrlen(descentFactors); i++)
	{
		//printf("descentFactor: ");fmpz_print(descentFactors[i]->p);printf("^%d\n",descentFactors[i]->exponent);		
	}
	return descentFactors;
}

void DescentFactor_UpdateFactors_FreeMem(fmpz_factor_t numeratorFactors, fmpz_factor_t denominatorFactors, DescentFactor *descentFactor, bool numOrden)
{
	fmpz_factor_struct *factorNum = numeratorFactors;
	fmpz_factor_struct *factorDen = denominatorFactors;

	if(numOrden == false)
	{
		//Working with DescentFactor den
		factorNum = denominatorFactors;
		factorDen = numeratorFactors;
	}
	for(size_t i = 0; i < arrlen(descentFactor); i++)
	{
		if(descentFactor[i]->exponent > 0)
		{
			//check if factor already exists
			bool found = false;
			for(slong j = 0; j < factorNum->num; j++)
			{
				if(fmpz_cmp(descentFactor[i]->p, &factorNum->p[j]) == 0)
				{
					//Increase exponent
					factorNum->exp[j] += abs(descentFactor[i]->exponent);
					found = true;
					break;
				}
			}
			if(found == false)
			{
				_fmpz_factor_append(factorNum, descentFactor[i]->p, abs(descentFactor[i]->exponent));			
			}
		}
		else
		{
			//check if factor already exists
			bool found = false;
			for(slong j = 0; j < factorDen->num; j++)
			{
				if(fmpz_cmp(descentFactor[i]->p, &factorDen->p[j]) == 0)
				{
					//Increase exponent
					factorDen->exp[j] += abs(descentFactor[i]->exponent);
					found = true;
					break;
				}
			}
			if(found == false)
			{
				_fmpz_factor_append(factorDen, descentFactor[i]->p, abs(descentFactor[i]->exponent));			
			}
		}
		DescentFactor_Destroy(descentFactor[i]);
	}
	arrfree(descentFactor);
}

void DescentSieve_FactorsToInt(fmpz_t result, fmpz_t prime, DescentFactor *descentFactor)
{
	fmpz_t resultNum,resultDen, temp0;
	fmpz_init(resultNum);fmpz_init(resultDen);fmpz_init(temp0);
	fmpz_set_ui(resultNum, 1);fmpz_set_ui(resultDen, 1);
	for(size_t i = 0; i < arrlen(descentFactor); i++)
	{
		fmpz_powm_ui(temp0, descentFactor[i]->p, abs(descentFactor[i]->exponent),prime);
		if(descentFactor[i]->exponent > 0)
		{
			fmpz_mul(resultNum, resultNum, temp0);
			fmpz_mod(resultNum, resultNum, prime);
		}
		else
		{
			fmpz_mul(resultDen, resultDen, temp0);
			fmpz_mod(resultDen, resultDen, prime);
		}
	}
	fmpz_invmod(resultDen,resultDen,prime);
	fmpz_mul(resultNum, resultNum, resultDen);	
	fmpz_mod(resultNum, resultNum, prime);
	//printf("DescentFactor res: ");fmpz_print(resultNum);printf("\n");
	fmpz_set(result, resultNum);
	fmpz_clear(resultNum);fmpz_clear(resultDen);fmpz_clear(temp0);
}

bool DescentSieve_ValidateReconstruction(fmpz_factor_t factorNum, fmpz_factor_t factorDen, fmpz_t target, fmpz_t prime)
{
	bool validReconstruction = false;
	fmpz_t resultNum,resultDen, temp0;
	fmpz_init(resultNum);fmpz_init(resultDen);fmpz_init(temp0);
	fmpz_set_ui(resultNum, 1);
	fmpz_set_ui(resultDen, 1);
	for(slong i = 0; i < factorNum->num; i++)
	{
		fmpz_powm_ui(temp0, &factorNum->p[i], factorNum->exp[i],prime);
		fmpz_mul(resultNum, resultNum, temp0);
		fmpz_mod(resultNum, resultNum, prime);
	}
	for(slong i = 0; i < factorDen->num; i++)
	{
		fmpz_powm_ui(temp0, &factorDen->p[i], factorDen->exp[i],prime);
		fmpz_mul(resultDen, resultDen, temp0);
		fmpz_mod(resultDen, resultDen, prime);
	}
	fmpz_invmod(resultDen,resultDen,prime);
	fmpz_mul(resultNum, resultNum, resultDen);	
	fmpz_mod(resultNum, resultNum, prime);
	validReconstruction = (fmpz_cmp(resultNum, target) == 0);
	//printf("Recon res: ");fmpz_print(resultNum);printf("\n");
	fmpz_clear(resultNum);fmpz_clear(resultDen);fmpz_clear(temp0);

	return validReconstruction;
}

void DescentSieveRun(fmpz_t target, fmpz_t prime, DescentSieve descentSieve)
{
	//Rational reconstruction
	uint64_t sieveX = 30000;//Max x value 300great
	uint64_t sieveY = 1;//100 great
	uint64_t maxYoffset = 10000;//100 great
	
	int numeratorSign = -1;
	int proved = 0;
	fmpz_t num, den, resultNum, resultDen, neg1;
	fmpz_init(num);fmpz_init(den);fmpz_init(resultNum);fmpz_init(resultDen);
	fmpz_init(neg1);fmpz_set_si(neg1, -1);
	fmpq_t rationalReconstruction;fmpq_init(rationalReconstruction);
	fmpz_factor_t factorNum;fmpz_factor_init(factorNum);
	fmpz_factor_t factorDen;fmpz_factor_init(factorDen);
	int reconstructionSuccess = fmpq_reconstruct_fmpz(rationalReconstruction, target, prime);
	if(reconstructionSuccess)
	{
		fmpz_abs(num, fmpq_numref(rationalReconstruction));
		fmpz_set(den, fmpq_denref(rationalReconstruction));
		//Find numerator sign
		if(fmpz_cmp_ui(fmpq_numref(rationalReconstruction), 0) < 0){numeratorSign = -1;}

		//Find Smooth factorization using descent sieve factorBaseBitLength
		int factorNumResult = fmpz_factor_smooth(factorNum, num, descentSieve->factorBaseBitLength, proved);
		int factorDenResult = fmpz_factor_smooth(factorDen, den, descentSieve->factorBaseBitLength, proved);
		
		//Sieve each factor bigger than our target length
		DescentFactor *descentNum = DescentSieve_Factors(factorNum, descentSieve,sieveX,sieveY,maxYoffset);
		DescentFactor *descentDen = DescentSieve_Factors(factorDen, descentSieve,sieveX,sieveY,maxYoffset);
		//DescentSieve_FactorsToInt(resultNum, prime, descentNum);
		//DescentSieve_FactorsToInt(resultDen, prime, descentDen);
		//Clear and init factorNum and factorDen
		fmpz_factor_clear(factorNum);fmpz_factor_clear(factorDen);
		fmpz_factor_init(factorNum);fmpz_factor_init(factorDen);
		if(numeratorSign == -1)
		{
			//Append neg 1
			_fmpz_factor_append(factorNum, neg1, 1);			
		}
		
		//Append DescentFactors and free DescentFactor memory
		DescentFactor_UpdateFactors_FreeMem(factorNum, factorDen, descentNum, true);		
		DescentFactor_UpdateFactors_FreeMem(factorNum, factorDen, descentDen, false);	
		
		//TestReconstruction
		bool validReconstruction = DescentSieve_ValidateReconstruction(factorNum, factorDen, target, prime);
		if(validReconstruction == false)
		{
			printf("validReconstruction is false for ");fmpz_print(target);printf("\n");
			assert(validReconstruction == true);
		}
		printf("Num: ");fmpz_print(fmpq_numref(rationalReconstruction));printf("\n");DescentSieve_PrintFactors(factorNum);
		printf("Den: ");fmpz_print(fmpq_denref(rationalReconstruction));printf("\n");DescentSieve_PrintFactors(factorDen);
	}
	fmpz_clear(num);fmpz_clear(den);fmpz_clear(neg1);
	fmpq_clear(rationalReconstruction);
	fmpz_factor_clear(factorNum);fmpz_factor_clear(factorDen);

}

void TestDescentSieve()
{
	fmpz_t target, prime;
	fmpz_init(target);fmpz_init(prime);
	char *targetStringBase10 = "314159265358979323846264338327950288419716939937510582097494459230781640628620899862";
	char *primeStringBase10  = "3108193808041961141219111205196826101966010119640309197118051941271219700607191207059";
	int factorBaseBitLength  = 25;

	fmpz_set_str(target, targetStringBase10,10);
	fmpz_set_str(prime, primeStringBase10,10);
	
	DescentSieve descentSieve= DescentSieve_Create(factorBaseBitLength,prime);
	DescentSieveRun(target, prime, descentSieve);
	
	DescentSieve_Destroy(descentSieve);
	fmpz_clear(target);fmpz_clear(prime);
}

int main()
{
	TestDescentSieve();
	flint_cleanup();
	return 0;
}
