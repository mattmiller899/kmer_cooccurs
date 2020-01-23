#include <stdio.h>
#include "uthash.h"
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <stdbool.h>
#include <libgen.h>
#include <string.h>
#include <stdlib.h>



#define BAD_NUC -128
#define BAD_AA -128
//#define kmer_size 3
//#define aa_size 1.0
//#define total_aas (int) pow(21.0, aa_size)

void test_cases();

int loop_counter = 0;
const bool exclude_plasmid = true;
int total_kmers;
int kmer_size, total_aas;
int aa_size;
struct kmer_struct *kmer2aa = NULL;
struct aa_count_struct *aa_counts = NULL;
int *local_aa_counts;
int *total_aa_counts;
char **kmer_ptrs;
const int read_size = 150;

char all_aas[21] = {'A', 'C', 'D', 'E', 'F', 'G',
                    'H', 'I', 'K', 'L', 'M', 'N',
                    'P', 'Q', 'R', 'S', 'T', 'V',
                    'W', 'Y', '*'};

char* codon_table[64] =
        {"K", "N", "K", "N",
         "T", "T", "T", "T",
         "R", "S", "R", "S",
         "I", "I", "M", "I",
         "Q", "H", "Q", "H",
         "P", "P", "P", "P",
         "R", "R", "R", "R",
         "L", "L", "L", "L",
         "E", "D", "E", "D",
         "A", "A", "A", "A",
         "G", "G", "G", "G",
         "V", "V", "V", "V",
         "*", "Y", "*", "Y",
         "S", "S", "S", "S",
         "*", "C", "W", "C",
         "L", "F", "L", "F"
        };

struct kmer_struct {
    char *kmer;            /* key */
    char *aa;
    UT_hash_handle hh; /* makes this structure hashable */
};
struct aa_count_struct {
    char *aa;            /* key */
    int *counts;
    UT_hash_handle hh; /* makes this structure hashable */
};
const char nucs[4] = "ACGT";
static const int8_t nuc_table[UCHAR_MAX+1] = {
        [0 ... 255] = BAD_NUC,  // ranges are a GNU extension
        ['A'] = 0,
        ['C'] = 1,
        ['G'] = 2,
        ['T'] = 3,
};

static const int8_t aa_table[UCHAR_MAX+1] = {
        [0 ... 255] = BAD_AA,
        ['A'] = 0,
        ['C'] = 1,
        ['D'] = 2,
        ['E'] = 3,
        ['F'] = 4,
        ['G'] = 5,
        ['H'] = 6,
        ['I'] = 7,
        ['K'] = 8,
        ['L'] = 9,
        ['M'] = 10,
        ['N'] = 11,
        ['P'] = 12,
        ['Q'] = 13,
        ['R'] = 14,
        ['S'] = 15,
        ['T'] = 16,
        ['V'] = 17,
        ['W'] = 18,
        ['Y'] = 19,
        ['*'] = 20,
};
static const int8_t codon_count_table[UCHAR_MAX+1] = {
        [0 ... 255] = BAD_AA,
        ['A'] = 4,
        ['C'] = 2,
        ['D'] = 2,
        ['E'] = 2,
        ['F'] = 2,
        ['G'] = 4,
        ['H'] = 2,
        ['I'] = 3,
        ['K'] = 2,
        ['L'] = 6,
        ['M'] = 1,
        ['N'] = 2,
        ['P'] = 4,
        ['Q'] = 2,
        ['R'] = 6,
        ['S'] = 6,
        ['T'] = 4,
        ['V'] = 4,
        ['W'] = 1,
        ['Y'] = 2,
        ['*'] = 3,
};
static const int8_t bp_table[UCHAR_MAX+1] = {
        [0 ... 255] = BAD_NUC,  // ranges are a GNU extension
        ['A'] = 'T',
        ['C'] = 'G',
        ['G'] = 'C',
        ['T'] = 'A',
};


// unsigned char* so high-ASCII -> 128..255, not negative,
// and works as an index into a 256 entry table
unsigned codon_to_idx_LUT(unsigned char *p) {
    unsigned idx = nuc_table[p[0]];
    idx = idx*4 + nuc_table[p[1]];
    idx = idx*4 + nuc_table[p[2]];
    return idx;
}

int ArgPos(char *str, int argc, char **argv) {
    int a;
    for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
            if (a == argc - 1) {
                printf("Argument missing for %s\n", str);
                exit(1);
            }
            return a;
        }
    return -1;
}


struct kmer_struct *find_kmer2aa(char *kmer_id) {
    struct kmer_struct *s;
    HASH_FIND_STR( kmer2aa, kmer_id, s );  /* s: output pointer */
    return s;
}
struct aa_count_struct *find_aa_count(char *aa_id) {
    struct aa_count_struct *a;
    HASH_FIND_STR( aa_counts, aa_id, a );  /* s: output pointer */
    //if(a == NULL) printf("%s not found\n", aa_id);
    return a;
}

//There will likely be memory leaks from mallocing a struct and then setting its parameters to another malloc'd variable
//But file scope errors make this the easiest solution
void init_kmer2aa(char *kmer_id, char *name) {
    if(find_kmer2aa(kmer_id) != NULL) return;
    struct kmer_struct *s;
    s = malloc(sizeof(struct kmer_struct));
    char *kmer = malloc(kmer_size+1);
    strcpy(kmer, kmer_id);
    char *aa = malloc(aa_size+1);
    strcpy(aa, name);
    s->kmer = kmer;
    s->aa = aa;
    //strcpy(s->kmer, kmer_id);
    //strcpy(s->aa, name);
    //printf("init struct kmer = %s aa = %s mem = %p\n", s->kmer, s->aa, &s);
    HASH_ADD_STR( kmer2aa, kmer, s );  /* id: name of key field */
}
void init_count(char *aa_id, int *tmp_counts) {
    if(find_aa_count(aa_id) != NULL) return;
    struct aa_count_struct *s;
    s = malloc(sizeof(struct aa_count_struct));
    char *aa = malloc(aa_size+1);
    int *counts = calloc(pow(21.0, aa_size), sizeof(int));
    strcpy(aa, aa_id);
    memcpy(counts, tmp_counts, sizeof(int)*21);
    s->aa = aa;
    s->counts = counts;
    //printf("init struct aa = %s count = %d mem = %p\n", s->aa, s->count, &s);
    HASH_ADD_STR( aa_counts, aa, s );  /* id: name of key field */
}
int *get_aa_counts(char *aa_id) {
    struct aa_count_struct *s = find_aa_count(aa_id);
    return s->counts;
}
void add_aa_counts(char *aa_id, int *add_counts) {
    struct aa_count_struct *s = find_aa_count(aa_id);
    for(int i = 0; i < total_aas; i++) {
        *(s->counts+i) += *(add_counts+i);
    }
}
char *convert_kmer2aa(char *kmer_id) {
    struct kmer_struct *s = find_kmer2aa(kmer_id);
    return s->aa;
}
unsigned aa_to_idx(char *tmp_aa) {
    unsigned total_idx = 0;
    int pow_exp = aa_size-1;
    for(int j = 0; j < aa_size; j++) {
        unsigned idx = aa_table[*(tmp_aa+j)];
        total_idx += (int) pow(21.0, pow_exp) * idx;
        pow_exp--;
    }
    return total_idx;
}

void idx_to_aa(int idx, char *tmp_aa) {
    int curr_pos, pow_exp, pow_val, dummy_count, tmp_idx;
    //memset(tmp_codon, 'X', kmer_size);
    curr_pos = 0;
    pow_exp = aa_size - 1;
    dummy_count = idx;
    //printf("starting count = %d\n", dummy_count);
    while (curr_pos < aa_size) {
        pow_val = (int) pow(21.0, (double) pow_exp);
        if (dummy_count < pow_val) tmp_aa[curr_pos++] = all_aas[0];
        else {
            tmp_idx = (int) (dummy_count / pow_val);
            tmp_aa[curr_pos++] = all_aas[tmp_idx];
            dummy_count -= pow_val * tmp_idx;
        }
        //printf("curr pos = %d pow_exp = %d dummy count = %d tmp codon = %s\n", curr_pos, pow_exp, dummy_count, tmp_codon);
        pow_exp--;
    }
}

void idx_to_kmer(int idx, char *tmp_codon) {
    int curr_pos, nuc_num, pow_exp, dummy_count, pow_val;
    curr_pos = 0;
    pow_exp = kmer_size - 1;
    dummy_count = idx;
    //printf("starting count = %d\n", dummy_count);
    while(curr_pos < kmer_size) {
        pow_val = (int) pow(4.0, (double)pow_exp);
        if(dummy_count < pow_val) tmp_codon[curr_pos++] = nucs[0];
        else {
            nuc_num = (int) (dummy_count / pow_val);
            tmp_codon[curr_pos++] = nucs[nuc_num];
            dummy_count -= pow_val * nuc_num;
        }
        //printf("curr pos = %d pow_exp = %d dummy count = %d tmp codon = %s\n", curr_pos, pow_exp, dummy_count, tmp_codon);
        pow_exp--;
    }
}

void generate_local_counts(char *aas, int start, int end, char *tmp_aa) {
    for(int i = start; i < end; i+=aa_size) {
        memcpy(tmp_aa, aas+i, aa_size);
        //printf("inside %d aa = %c\n", i, tmp_aa);
        //Calculate index for aa_counts
        unsigned idx = aa_to_idx(tmp_aa);
        //printf("idx = %d\n", idx);
        //if(idx == 119) printf("YYYYYEEEEEAHHHHH %d mem = %p\n", i, &aas+i);
        //printf("total_idx = %d aa = %s\n", total_idx, tmp_aa);
        *(local_aa_counts + idx) += 1;
    }
}
void check_local_counts(int old) {
    int sum = 0;
    bool error_found = false;
    for(int i = 0; i < total_aas; i++) {
        if(*(local_aa_counts+i) < 0) {
            printf("ERROR: old %d local_aa_counts contains value less than 0: idx %d = %d\n", old, i, *(local_aa_counts+i));
            error_found = true;
        }
        sum += *(local_aa_counts+i);
    }
    if(sum != 2*read_size/aa_size + 1) {
        printf("ERROR: sum of local counts != (2*read size)/aa_size: %d != %d\n", 2*read_size/aa_size+1, sum);
        error_found = true;
    }
    if(error_found) exit(1);
}
void generate_cooccurs(char *aas, int num_aas) {
    int beg = num_aas - read_size;
    int end = read_size;
    memset(local_aa_counts, 0, sizeof(int) * total_aas);
    printf("in cooccur num kmers = %d\n", num_aas);
    /*
    char *first_hun = calloc(100, sizeof(char));
    strncpy(first_hun, aas, 100);
    printf("first hundred = %s\n", first_hun);
    */
    //Generate initial local aa counts
    char *tmp_aa = malloc(aa_size+1);
    tmp_aa[(int) aa_size] = '\0';
    generate_local_counts(aas, 0, end+1, tmp_aa);
    generate_local_counts(aas, beg, num_aas, tmp_aa);
    int sum = 0;
    for(int i = 0; i < total_aas; i++) {
        sum += *(local_aa_counts+i);
        //if(*(local_aa_counts+i) < 0) printf("ERROR: local_aa_counts contains value less than 0: idx %d = %d\n", i, *(local_aa_counts+i));
    }
    if(sum != 2*read_size/aa_size + 1) printf("ERROR: sum of local counts != (2*read size)/aa_size + 1: %d != %d\n", 2*read_size/aa_size + 1, sum);
    //Scan over the AAs, adding the local AAs to the AA hashmap and moving curr, beg, and end
    for(int i = 0; i < num_aas - aa_size + 1; i+=aa_size) {
        //printf("tmp_aa = %s\n", tmp_aa);
        //Add current local counts to the current A
        //Subtract AA at "i" as it doesnt cooccur with itself, then readd back
        memcpy(tmp_aa, aas+i, aa_size);
        unsigned aa_idx = aa_to_idx(tmp_aa);
        local_aa_counts[aa_idx]--;
        add_aa_counts(tmp_aa, local_aa_counts);
        local_aa_counts[aa_idx]++;
        //printf("added\n");
        //Shift beg and end, adding and subtracting AAs from local
        memcpy(tmp_aa, aas+beg, aa_size);
        aa_idx = aa_to_idx(tmp_aa);
        //if(aa_idx == 119) printf("SUBBING %d mem = %p\n", beg, &aas+beg);
        //if(aa_idx == 5 && loop_counter == 2) printf("subtracting 1 count = %d\n", *(local_aa_counts+5) - 1);
        local_aa_counts[aa_idx]--;
        //Subtract AA size (so that if, for example, aa_size = 2 and end is num_aas - 1, the null character isn't in tmp_aa
        //Add +1 because
        beg = (int)(beg + aa_size) % (int) (num_aas);
        end = (int)(end + aa_size) % (int) (num_aas);
        if(end == num_aas - 1) printf("DANGER\n");
        memcpy(tmp_aa, aas+end, aa_size);
        aa_idx = aa_to_idx(tmp_aa);
        //if(aa_idx == 119) printf("ADDING %d mem = %p\n", end, &aas+end);
        //if(aa_idx == 5 && loop_counter == 2) printf("adding 1 count = %d\n", *(local_aa_counts+5) + 1);
        local_aa_counts[aa_idx]++;
        //TODO COMMENT OUT
        check_local_counts(i);
        //printf("shifted\n");
        /*
        int *tmp = get_aa_counts(tmp_aa);
        for(int j = 0; j < 21; j++) {
            printf("curr aa = %c count %d = %d\n", *(aas+i), j, *(tmp+j));
        }
        */
    }
}

void write_cooccurs(char *out_fp) {
    FILE *out = fopen(out_fp, "w+");
    //Print header
    char *tmp_aa = malloc(aa_size+1);
    tmp_aa[(int)aa_size] = '\0';
    fprintf(out, "AAs");
    int total_count = 0;
    while(total_count < total_aas) {
        fprintf(out, ",");
        //memset(tmp_codon, 'X', kmer_size);
        idx_to_aa(total_count, tmp_aa);
        total_count++;
        fprintf(out, "%s", tmp_aa);
    }
    fprintf(out, "\n");
    /*
    for(int i = 0; i < total_aas - 1; i++) {
        //memcpy(tmp_aa, )
        fprintf(out, "%c,", all_aas[i]);
    }
    fprintf(out, "%c\n", all_aas[20]);
    */
    //Print AA then its cooccurs
    for(int i = 0; i < total_aas; i++) {
        //TODO PROLLY WONT WORK FOR 6mer
        idx_to_aa(i, tmp_aa);
        struct aa_count_struct *a = find_aa_count(tmp_aa);
        fprintf(out, "%s", tmp_aa);
        for(int j = 0; j < total_aas; j++) {
            fprintf(out, ",%d", a->counts[j]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
    free(tmp_aa);
}

void norm_aa_counts() {
    char *tmp_aa = malloc(aa_size+1);
    tmp_aa[(int)aa_size] = '\0';
    char *tmp_aa2 = malloc(aa_size+1);
    tmp_aa2[(int)aa_size] = '\0';
    for(int i = 0; i < total_aas; i++) {
        //TODO PROLLY WONT WORK FOR 6mer
        idx_to_aa(i, tmp_aa);
        int *tmp_counts = get_aa_counts(tmp_aa);
        for (int j = 0; j < total_aas; j++) {
            idx_to_aa(i, tmp_aa2);
            int num_codons = 1;
            //Normalizing by # of codons for 1st AA times # of codons for 2nd AA
            for(int k = 0; k < aa_size; k++) {
                num_codons = num_codons * codon_count_table[*(tmp_aa+k)] * codon_count_table[*(tmp_aa2+k)];
            }
            //TODO rounding "errors"
            tmp_counts[j] = (int) tmp_counts[j] / num_codons;
        }
    }
    free(tmp_aa);
    free(tmp_aa2);
}

//From https://stackoverflow.com/questions/2736753/how-to-remove-extension-from-file-name
char *remove_ext (char* myStr, char extSep, char pathSep) {
    char *retStr, *lastExt, *lastPath;
    if (myStr == NULL) return NULL;
    if ((retStr = malloc (strlen (myStr) + 1)) == NULL) return NULL;
    strcpy (retStr, myStr);
    lastExt = strrchr (retStr, extSep);
    lastPath = (pathSep == 0) ? NULL : strrchr (retStr, pathSep);
    if (lastExt != NULL) {
        if (lastPath != NULL) {
            if (lastPath < lastExt) {
                *lastExt = '\0';
            }
        } else {
            *lastExt = '\0';
        }
    }
    return retStr;
}

char *my_itoa(int num, char *str)
{
    if(str == NULL)
    {
        return NULL;
    }
    sprintf(str, "%d", num);
    return str;
}


int main(int argc, char **argv) {
    printf("Hello, World!\n");
    char * line = calloc(10000000, 1);
    size_t len = 0;
    ssize_t chars_read;
    FILE *fin;
    //Initialize file TODO GET FROM ARG
    char fp[500];
    char out_dir[500];
    int i;
    if ((i = ArgPos((char *)"-in_file", argc, argv)) > 0) strcpy(fp, argv[i + 1]);
    if ((i = ArgPos((char *)"-out_dir", argc, argv)) > 0) strcpy(out_dir, argv[i + 1]);
    if ((i = ArgPos((char *)"-kmer", argc, argv)) > 0) kmer_size = atoi(argv[i + 1]);
    aa_size = kmer_size / 3;
    total_aas = (int) pow(21.0, (double) aa_size);
    total_kmers = (int) pow(4.0, (double) kmer_size);
    local_aa_counts = (int*) calloc(total_aas, sizeof(int));
    char *fn = remove_ext(basename(fp), '.', '/');
    char *raw_fp = calloc(500, 1);
    char *norm_fp = calloc(500, 1);
    char *kmer_str = calloc(2, 1);
    my_itoa(kmer_size, kmer_str);
    strcat(kmer_str, "mer.csv");
    printf("kmerstr = %s\n", kmer_str);
    strcpy(raw_fp, out_dir);
    strcat(raw_fp, "/");
    strcat(raw_fp, fn);
    strcat(raw_fp, "_cooccurs_raw_");
    strcat(raw_fp, kmer_str);
    printf("rawfp = %s\n", raw_fp);
    strcpy(norm_fp, out_dir);
    strcat(norm_fp, "/");
    strcat(norm_fp, fn);
    strcat(norm_fp, "_cooccurs_norm_");
    strcat(norm_fp, kmer_str);
    printf("norm = %s\n", norm_fp);


    //Initialize hashtables
    /*
    if(aa_size == 1.0) {
        char *tmp_codon = malloc(kmer_size+1);
        int *tmp_counts = calloc(21, sizeof(int));
        for(int i = 0; i < 4; i++) {
            tmp_codon[0] = nucs[i];
            for(int j = 0; j < 4; j++) {
                tmp_codon[1] = nucs[j];
                for(int k = 0; k < 4; k++) {
                    tmp_codon[2] = nucs[k];
                    unsigned tmp_idx = codon_to_idx_LUT(tmp_codon);
                    char *tmp_aa = codon_table[tmp_idx];
                    //printf("codon = %s aa = %s\n", tmp_codon, tmp_aa);
                    init_kmer2aa(tmp_codon, tmp_aa);
                    // TODO CHECK MEMORY FOR 3mers
                    struct kmer_struct *s = find_kmer2aa(tmp_codon);
                    printf("kmer2aa kmer = %s aa = %s mem = %p\n", s->kmer, s->aa, &s);
                    init_count(tmp_aa, tmp_counts);
                    struct aa_count_struct *a = find_aa_count(codon_table[tmp_idx]);
                    printf("aa_counts aa = %s count = %d mem = %p\n", a->aa, a->counts[0], &a);
                }
            }
        }
        free(tmp_codon);
        free(tmp_counts);
    }
    */
    //else {
    int total_count = 0;
    int *tmp_counts = calloc(total_aas, sizeof(int));
    char *tmp_codon = malloc(kmer_size+1);
    *(tmp_codon+kmer_size) = '\0';
    char *tmp_aa = calloc(aa_size+1, 1);
    tmp_aa[(int) aa_size] = '\0';
    char *one_codon = calloc(4, 1);
    *(one_codon + 3) = '\0';
    unsigned tmp_idx;
    while(total_count < total_kmers) {
        //memset(tmp_codon, 'X', kmer_size);
        idx_to_kmer(total_count, tmp_codon);
        total_count++;
        //printf("%s\n", tmp_codon);
        for(int i = 0; i < kmer_size; i+= 3) {
            memcpy(one_codon, tmp_codon + i, 3);
            tmp_idx = codon_to_idx_LUT(one_codon);
            *(tmp_aa + (i / 3)) = *codon_table[tmp_idx];
        }
        //printf("aa = %s\n", tmp_aa);
        init_kmer2aa(tmp_codon, tmp_aa);
        init_count(tmp_aa, tmp_counts);
        /*
        struct kmer_struct *s = find_kmer2aa(tmp_codon);
        printf("kmer2aa kmer = %s aa = %s mem = %p\n", s->kmer, s->aa, &s);
        struct aa_count_struct *a = find_aa_count(tmp_aa);
        printf("aa_counts aa = %s count = %d mem = %p\n", a->aa, a->counts[0], &a);
        struct aa_count_struct *b = find_aa_count(prev_aa);
        printf("prev aa_counts aa = %s count = %d mem = %p\n", b->aa, b->counts[0], &b);
        memcpy(prev_aa, tmp_aa, 3);
        */
    }
    free(tmp_counts);
    //}
    //Read the fasta header
    fin = fopen(fp, "rb");
    chars_read = getline(&line, &len, fin);
    //kmer_ptrs = calloc(kmer_size, sizeof(char*));
    bool skip = false;
    //Read lines, skipping headers (and sequence if its from a plasmid), generating AAs and finding cooccurs
    while ((chars_read = getline(&line, &len, fin)) != -1) {
        if(line[0] == '>') {
            printf("header found %s\n", line);
            skip = false;
            if(strstr(line, "plasmid"))
            {
                printf("header with plasmid found %s\n", line);
                if(exclude_plasmid) {
                    skip = true;
                }
            }
            continue;
        }
        if(skip) continue;
        int max_num_kmers = ceil(chars_read / 3);
        for(int i = 0; i < kmer_size; i++) {
            loop_counter = i;
            char *kmer_ptr = calloc(max_num_kmers, sizeof(char));
            int counted_num_kmers = 0;
            printf("here %d\n", i);
            //memcpy(kmer_ptrs+(i*sizeof(char*)), &kmer_ptr, sizeof(char*));
            for(int j = i; j < chars_read - kmer_size; j+=kmer_size) {
                memcpy(tmp_codon, line+j, kmer_size);
                //printf("tmp codon = %s\n", tmp_codon);
                tmp_aa = convert_kmer2aa(tmp_codon);
                //printf("tmp aa = %s\n", tmp_aa);
                memcpy(kmer_ptr+counted_num_kmers, tmp_aa, aa_size);
                counted_num_kmers+=aa_size;
                //printf("kmer_ptr = %s\n", kmer_ptr);
            }
            //printf("done with %d\n", i);
            //printf("kmer_ptrs %d = %s\n", i, kmer_ptrs[i]);
            //printf("kmer_ptr %d = %s\n", i, kmer_ptr);
            printf("done %d\n", i);
            generate_cooccurs(kmer_ptr, counted_num_kmers);
            printf("finished cooccurs %d\n", i);
            free(kmer_ptr);
        }
        printf("finished with ORFS\n");
    }
    printf("Done calculating cooccurs, writing and normalizing\n");
    fclose(fin);
    free(tmp_codon);
    //Write the cooccur stats to two file, 1 raw 1 normed on # of codons/AA (W has 1 codon, L has 6 codons)
    write_cooccurs(raw_fp);
    norm_aa_counts();
    write_cooccurs(norm_fp);

}
/*
void test_cases() {
    //Test 1: aa table matches for 3mers
    char tmp_codon[3] = "AAA";
    if(strcmp(codon_table[codon_to_idx_LUT(tmp_codon)], "K") != 0) printf("ERROR IN TEST CASE 1: table %s != tested %s\n",
                                                                       codon_table[codon_to_idx_LUT(tmp_codon)], "K");
    memcpy(tmp_codon, "CAG", 3);
    if(strcmp(codon_table[codon_to_idx_LUT(tmp_codon)], "Q") != 0) printf("ERROR IN TEST CASE 1: table %s != tested %s\n",
                                                                       codon_table[codon_to_idx_LUT(tmp_codon)], "Q");

    //Test 2: aa_counts increment correctly
    add_aa_counts(codon_table[4], 2);
    if(get_aa_counts(codon_table[4]) != 2) printf("ERROR IN TEST CASE 2: add_aa_count failed %d != %d", codon_table[4], 2);

    struct kmer_struct *s = find_kmer2aa("AAA");
    printf("struct id = %s\n", s->kmer);
    s = find_kmer2aa("CAG");
    printf("struct id = %s aa = %s\n", s->kmer, s->aa);
}
*/
