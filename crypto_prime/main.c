#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#define MAX_PRIME 100000000

typedef union prime_data
{
    mpz_t mpz;
    unsigned long simple;
}prime_data_u;

typedef struct prime_numbers
{
    prime_data_u* array;
    unsigned int prime_numbers_nr;
}prime_numbers_t;
static prime_numbers_t * prime_number_array;

typedef struct bezout
{
    mpz_t u;
    mpz_t v;
} bezout_t;

typedef struct chunk
{
    prime_numbers_t * prime_number_array;
    struct chunk* next;
} chunk_t;

chunk_t* head;
chunk_t* current;


struct primes
{
  unsigned int prime;
  int rem;
};

struct primes *primes;
unsigned long n_primes;

static mpz_t settings_from;
static mpz_t settings_to;

void find_primes (unsigned char *, mpz_t, unsigned long, mpz_t);
void sieve_region (unsigned char *, mpz_t, unsigned long);
void make_primelist (unsigned long);

int flag_print = 0;
int flag_count = 0;
int flag_maxgap = 0;
unsigned long maxgap = 0;
unsigned long total_primes = 0;

void init_prime_array(prime_numbers_t** prime_numbers)
{

*prime_numbers = malloc(sizeof(prime_numbers_t));
(*prime_numbers)->array = malloc(MAX_PRIME * sizeof(*(*prime_numbers)->array));
if ((*prime_numbers)->array == NULL)
{
    printf("malloc pb\n");    
    exit(1);
}
prime_number_array= *prime_numbers;
}

void
report (mpz_t prime)
{

  total_primes += 1;
 if (prime_number_array->prime_numbers_nr == MAX_PRIME)
  {
      printf("new chunk\n");
    chunk_t*  new_chunk = malloc(sizeof(chunk_t));
    current->next = new_chunk;
    init_prime_array(&new_chunk->prime_number_array);
    new_chunk->next = NULL;

  }
    if (mpz_fits_ulong_p(prime))
    {
        prime_number_array->array[prime_number_array->prime_numbers_nr].simple = mpz_get_ui(prime);
    }
    else
    {
        printf("need mpz\n");
        prime_number_array->array[prime_number_array->prime_numbers_nr].simple = 0;
        mpz_init(prime_number_array->array[prime_number_array->prime_numbers_nr].mpz);
        mpz_set(prime_number_array->array[prime_number_array->prime_numbers_nr].mpz, prime);
    }
    prime_number_array->prime_numbers_nr++;
  if (flag_print)
    {
      mpz_out_str (stdout, 10, prime);
      printf ("\n");
    }

}

void
generate_primes ()
{
  mpz_t fr, to;
  mpz_t fr2, to2;
  unsigned long sieve_lim;
  unsigned long est_n_primes;
  unsigned char *s;
  mpz_t tmp;
  mpz_t siev_sqr_lim;

  mpz_init (fr);
  mpz_init (to);
  mpz_init (fr2);
  mpz_init (to2);

      mpz_set (fr, settings_from);
      mpz_set (to, settings_to);

  mpz_set (fr2, fr);
  if (mpz_cmp_ui (fr2, 3) < 0)
    {
      mpz_set_ui (fr2, 2);
      report (fr2);
      mpz_set_ui (fr2, 3);
    }
  mpz_setbit (fr2, 0);                          /* make odd */
  mpz_sub_ui (to2, to, 1);
  mpz_setbit (to2, 0);                          /* make odd */

  mpz_init (tmp);
  mpz_init (siev_sqr_lim);

  mpz_sqrt (tmp, to2);
#define SIEVE_LIMIT 10000000
  if (mpz_cmp_ui (tmp, SIEVE_LIMIT) < 0)
    {
      sieve_lim = mpz_get_ui (tmp);
    }
  else
    {
      sieve_lim = SIEVE_LIMIT;
      mpz_sub (tmp, to2, fr2);
      if (mpz_cmp_ui (tmp, sieve_lim) < 0)
        sieve_lim = mpz_get_ui (tmp);   /* limit sieving for small ranges */
    }
  mpz_set_ui (siev_sqr_lim, sieve_lim + 1);
  mpz_mul_ui (siev_sqr_lim, siev_sqr_lim, sieve_lim + 1);

  est_n_primes = (size_t) (sieve_lim / log((double) sieve_lim) * 1.13) + 10;
  primes = malloc (est_n_primes * sizeof primes[0]);
  make_primelist (sieve_lim);
  assert (est_n_primes >= n_primes);

#if DEBUG
  printf ("sieve_lim = %lu\n", sieve_lim);
  printf ("n_primes = %lu (3..%u)\n",
          n_primes, primes[n_primes - 1].prime);
#endif

#define S (1 << 15)             /* FIXME: Figure out L1 cache size */
  s = malloc (S/2);
  while (mpz_cmp (fr2, to2) <= 0)
    {
      unsigned long rsize;
      rsize = S;
      mpz_add_ui (tmp, fr2, rsize);
      if (mpz_cmp (tmp, to2) > 0)
        {
          mpz_sub (tmp, to2, fr2);
          rsize = mpz_get_ui (tmp) + 2;
        }
#if DEBUG
      printf ("Sieving region ["); mpz_out_str (stdout, 10, fr2);
      printf (","); mpz_add_ui (tmp, fr2, rsize - 2);
      mpz_out_str (stdout, 10, tmp); printf ("]\n");
#endif
      sieve_region (s, fr2, rsize);
      find_primes (s, fr2, rsize / 2, siev_sqr_lim);

      mpz_add_ui (fr2, fr2, S);
    }
  free (s);
    free(primes);

}

/* Find primes in region [fr,fr+rsize).  Requires that fr is odd and that
   rsize is even.  The sieving array s should be aligned for "long int" and
   have rsize/2 entries, rounded up to the nearest multiple of "long int".  */
void
sieve_region (unsigned char *s, mpz_t fr, unsigned long rsize)
{
  unsigned long ssize = rsize / 2;
  unsigned long start, start2, prime;
  unsigned long i;
  mpz_t tmp;

  mpz_init (tmp);

#if 0
  /* initialize sieving array */
  for (ii = 0; ii < (ssize + sizeof (long) - 1) / sizeof (long); ii++)
    ((long *) s) [ii] = ~0L;
#else
  {
    long k;
    long *se = (long *) (s + ((ssize + sizeof (long) - 1) & -sizeof (long)));
    for (k = -((ssize + sizeof (long) - 1) / sizeof (long)); k < 0; k++)
      se[k] = ~0L;
  }
#endif

  for (i = 0; i < n_primes; i++)
    {
      prime = primes[i].prime;

      if (primes[i].rem >= 0)
        {
          start2 = primes[i].rem;
        }
      else
        {
          mpz_set_ui (tmp, prime);
          mpz_mul_ui (tmp, tmp, prime);
          if (mpz_cmp (fr, tmp) <= 0)
            {
              mpz_sub (tmp, tmp, fr);
              if (mpz_cmp_ui (tmp, 2 * ssize) > 0)
                break;          /* avoid overflow at next line, also speedup */
              start = mpz_get_ui (tmp);
            }
          else
            {
              start = (prime - mpz_tdiv_ui (fr, prime)) % prime;
              if (start % 2 != 0)
                start += prime;         /* adjust if even divisible */
            }
          start2 = start / 2;
        }

#if 0
      for (ii = start2; ii < ssize; ii += prime)
        s[ii] = 0;
      primes[i].rem = ii - ssize;
#else
      {
        long k;
        unsigned char *se = s + ssize; /* point just beyond sieving range */
        for (k = start2 - ssize; k < 0; k += prime)
          se[k] = 0;
        primes[i].rem = k;
      }
#endif
    }
  mpz_clear (tmp);
}

/* Find primes in region [fr,fr+rsize), using the previously sieved s[].  */
void
find_primes (unsigned char *s, mpz_t  fr, unsigned long ssize,
             mpz_t siev_sqr_lim)
{
  unsigned long j, ij;
  mpz_t tmp;

  mpz_init (tmp);
  for (j = 0; j < (ssize + sizeof (long) - 1) / sizeof (long); j++)
    {
      if (((long *) s) [j] != 0)
        {
          for (ij = 0; ij < sizeof (long); ij++)
            {
              if (s[j * sizeof (long) + ij] != 0)
                {
                  if (j * sizeof (long) + ij >= ssize)
                    goto out;
                  mpz_add_ui (tmp, fr, (j * sizeof (long) + ij) * 2);
                  if (mpz_cmp (tmp, siev_sqr_lim) < 0 ||
                      mpz_probab_prime_p (tmp, 10))
                    report (tmp);
                }
            }
        }
    }
 out:
  mpz_clear (tmp);
}

/* Generate a list of primes and store in the global array primes[].  */
void
make_primelist (unsigned long maxprime)
{
#if 1
  unsigned char *s;
  unsigned long ssize = maxprime / 2;
  unsigned long i, ii, j;

  s = malloc (ssize);
  memset (s, ~0, ssize);
  for (i = 3; ; i += 2)
    {
      unsigned long isqr = i * i;
      if (isqr >= maxprime)
        break;
      if (s[i * i / 2 - 1] == 0)
        continue;                               /* only sieve with primes */
      for (ii = i * i / 2 - 1; ii < ssize; ii += i)
        s[ii] = 0;
    }
  n_primes = 0;
  for (j = 0; j < ssize; j++)
    {
      if (s[j] != 0)
        {
          primes[n_primes].prime = j * 2 + 3;
          primes[n_primes].rem = -1;
          n_primes++;
        }
    }
  /* FIXME: This should not be needed if fencepost errors were fixed... */
  if (primes[n_primes - 1].prime > maxprime)
    n_primes--;
  free (s);
#else
  unsigned long i;
  n_primes = 0;
  for (i = 3; i <= maxprime; i += 2)
    {
      if (i < 7 || (i % 3 != 0 && i % 5 != 0 && i % 7 != 0))
        {
          primes[n_primes].prime = i;
          primes[n_primes].rem = -1;
          n_primes++;
        }
    }
#endif
}

void serialize_prime_numbers(const char*  filename )
{
    unsigned int i;
    FILE* fp = fopen(filename, "w");
    for (i = 0; i < prime_number_array->prime_numbers_nr; i++)
    {
        if (prime_number_array->array[i].simple ==0 )
        {
            mpz_out_str (fp, 10, prime_number_array->array[i].mpz);
        }
        else
            fprintf(fp, "%lu", prime_number_array->array[i].simple);            

        fprintf(fp, "\n"); 
    }
    fclose(fp);
}

void deserialize_prime_numbers(const char*  filename)
{
    FILE* fp = fopen(filename, "r");
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    prime_number_array->prime_numbers_nr = 0;
    mpz_t tmp;
    mpz_init(tmp);
    while ((read = getline(&line, &len, fp)) != -1) {

        mpz_set_str(tmp, line, 10);
    if (mpz_fits_ulong_p(tmp))
    {
        prime_number_array->array[prime_number_array->prime_numbers_nr].simple = mpz_get_ui(tmp);
    }
    else
    {
        printf("need mpz\n");
        prime_number_array->array[prime_number_array->prime_numbers_nr].simple = 0;
        mpz_init(prime_number_array->array[prime_number_array->prime_numbers_nr].mpz);
        mpz_set(prime_number_array->array[prime_number_array->prime_numbers_nr].mpz, tmp);
    }
        prime_number_array->prime_numbers_nr++;


    }
    fclose(fp);
}
/*
void filter_prime_number(unsigned int digit_nr, prime_numbers_t* prime_numbers_array, prime_numbers_t* res)
{
    unsigned int i;
    res->prime_numbers_nr = 0;
    unsigned int boundary = 0;
    if (digit_nr > 0)
    {
        boundary = (int) pow((double)10, (double)digit_nr - 1);

    }
    for (i = 0; i <  prime_numbers_array->prime_numbers_nr; i++)
    {
        if (prime_numbers_array->array[i] >= boundary)
        {
            res->array[res->prime_numbers_nr++] = prime_numbers_array->array[i];
        }
    }
}
*/
void bezout(mpz_t  a, mpz_t b, bezout_t* res)
{
    mpz_t r,rp,u,v,up,vp;
    mpz_init(r);
    mpz_init(rp);
    mpz_init(u);
    mpz_init(v);
    mpz_init(up);
    mpz_init(vp);
    mpz_set(r,a);

    mpz_set(rp,b);
    mpz_set_ui(u,1);
    mpz_set_ui(v,0);
    mpz_set_ui(up,0);
    mpz_set_ui(vp,1);

    mpz_t tmp, tmp2;
    mpz_init(tmp);
    mpz_init(tmp2);
    while (mpz_cmp_ui(rp, 0) != 0)
    {
        mpz_t q,rs,us,vs;

    mpz_init(q);
    mpz_init(rs);
    mpz_init(us);
    mpz_init(vs);
    mpz_mod(tmp,r, rp);
    mpz_fdiv_q(q, r, rp);
    mpz_set(rs,r);
    mpz_set(us,u);
    mpz_set(vs,v);
    mpz_set(r,rp);
    mpz_set(u,up);
    mpz_set(v,vp);
    mpz_mul(tmp,q, rp);
    mpz_sub(rp, rs, tmp);

    mpz_mul(tmp,q, up);
    mpz_sub(up, us, tmp);

    mpz_mul(tmp,q, vp);
    mpz_sub(vp, vs, tmp);
    mpz_clear (q);
    mpz_clear (rs);
    mpz_clear (us);
    mpz_clear (vs);

    }
    mpz_init(res->u);
    mpz_init(res->v);
    mpz_set(res->u,u);
    mpz_set(res->v,v);
}


void crack(const char * nb)
{
    mpz_t n, tmp;
    mpz_init(n);
    mpz_init(tmp);
    mpz_set_str(n, nb, 10);
    unsigned int i;
        mpz_t res_mod;
        mpz_init(res_mod);
    current = head;
    while (current)
    {
        prime_number_array = current->prime_number_array;
    for (i = 0; i < prime_number_array->prime_numbers_nr; i++)
    {
        if (prime_number_array->array[i].simple)
        {
            mpz_set_ui (tmp, prime_number_array->array[i].simple);
        }
        else
        {
            mpz_set(tmp, prime_number_array->array[i].mpz);
        }
        mpz_mod (res_mod,n,tmp);
        if (mpz_cmp_ui(res_mod, 0) == 0)
        {
            mpz_t p;
           mpz_t q;
            mpz_init (q);
            mpz_init (p);
           mpz_set(p,tmp);
           mpz_divexact(q,n,p);
            printf("p = ");
            mpz_out_str (stdout, 10, p);
            printf ("\n");
            printf("q = ");
            mpz_out_str (stdout, 10, q);
            printf ("\n");
            bezout_t res_bezout;
            mpz_t p1, q1, phi;
            mpz_init (q1);
            mpz_init (p1);
            mpz_init (phi); 
            mpz_sub_ui (p1, p, 1);
            mpz_sub_ui (q1, q, 1);
            mpz_mul(phi, p1, q1);
            bezout(n, phi, &res_bezout);
            printf("u = ");
            mpz_out_str (stdout, 10, res_bezout.u);
            printf ("\n");
            printf("v = ");
            mpz_out_str (stdout, 10, res_bezout.v);
            printf ("\n");

            return;
        }
    }
    current = current->next;
    }
}




int main(int argc, char *argv[])
{
clock_t begin, end;
double time_spent;
head = current =  malloc(sizeof(chunk_t));
init_prime_array(&head->prime_number_array);
head->next = NULL;

mpz_init(settings_from);
mpz_init(settings_to);
char from_s[2048];
char to_s[2048];
char n_s[2048];
printf("settings : bornes inf pour primes : \n");
scanf("%s", from_s);


printf("settings : bornes sup pour primes : \n");
scanf("%s", to_s);


    mpz_set_str(settings_from, from_s, 10);

    mpz_set_str(settings_to, to_s, 10);
begin = clock();
generate_primes();
end = clock();
time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Find primes from %s to %s : %f seconds\n", from_s, to_s, time_spent);

//deserialize_prime_numbers("out.prime");
while(1)
{   
    printf("saisir n : \n");
   scanf("%s", n_s);
   begin = clock(); 
    crack(n_s);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time : %f seconds\n", time_spent);
}
    /*deserialize_prime_numbers("out.prime", prime_numbers_array);
    prime_numbers_t* prime_numbers_array_filtered =  malloc (sizeof(prime_numbers_t));
    filter_prime_number(0, prime_numbers_array, prime_numbers_array_filtered);
    unsigned long long n = strtoull(argv[1], NULL, 0);
    crack(n, prime_numbers_array_filtered);*/
    return 0;
}
