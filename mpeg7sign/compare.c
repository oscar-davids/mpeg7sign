#include "signature.h"

#define HOUGH_MAX_OFFSET 90
#define MAX_FRAMERATE 60

#define DIR_PREV 0
#define DIR_NEXT 1
#define DIR_PREV_END 2
#define DIR_NEXT_END 3

#define STATUS_NULL 0
#define STATUS_END_REACHED 1
#define STATUS_BEGIN_REACHED 2


static void fill_l1distlut(uint8_t lut[])
{
	for (int i = 0, count = 0; i < 242; ++i) {
		for (int j = i + 1; j < 243; ++j, ++count) {
			/* ternary distance between i and j */
			uint8_t dist = 0;
			int tmp_i = i, tmp_j = j;
			do {
				dist += FFABS((tmp_j % 3) - (tmp_i % 3));
				tmp_j /= 3;
				tmp_i /= 3;
			} while (tmp_i > 0 || tmp_j > 0);
			lut[count] = dist;
		}
	}
}

static unsigned int intersection_word(const uint8_t *first,	const uint8_t *second)
{
	unsigned int val = 0;
	for (unsigned int i = 0; i < 28; i += 4) {
		val += av_popcount((first[i] & second[i]) << 24 |
			(first[i + 1] & second[i + 1]) << 16 |
			(first[i + 2] & second[i + 2]) << 8 |
			(first[i + 3] & second[i + 3]));
	}
	val += av_popcount((first[28] & second[28]) << 16 |
		(first[29] & second[29]) << 8 |
		(first[30] & second[30]));
	return val;
}

static unsigned int union_word(const uint8_t *first, const uint8_t *second)
{
	unsigned int val = 0;
	for (unsigned int i = 0; i < 28; i += 4) {
		val += av_popcount((first[i] | second[i]) << 24 |
			(first[i + 1] | second[i + 1]) << 16 |
			(first[i + 2] | second[i + 2]) << 8 |
			(first[i + 3] | second[i + 3]));
	}
	val += av_popcount((first[28] | second[28]) << 16 |
		(first[29] | second[29]) << 8 |
		(first[30] | second[30]));
	return val;
}


//__attribute__((optimize("unroll-loops")))
static unsigned int get_l1dist(SignatureContext *sc, const uint8_t *first, const uint8_t *second)
{
	unsigned int dist = 0;

	for (unsigned int i = 0; i < SIGELEM_SIZE / 5; ++i) {
		// original code kept for clarity purpose
		uint8_t f = first[i];
		uint8_t s = second[i];
		if (f != s) {
			if (f > s)
				// little variation of gauss sum formula
				dist += sc->l1distlut[243 * 242 / 2 - (243 - s)*(242 - s) / 2 + f - s - 1];
			else
				dist += sc->l1distlut[243 * 242 / 2 - (243 - f)*(242 - f) / 2 + s - f - 1];
		}
	}

	return dist;
}


static int get_jaccarddist(SignatureContext *sc, CoarseSignature *first, CoarseSignature *second)
{
	int composdist = 0, cwthcount = 0;
	for (int i = 0; i < 5; i++) {
		unsigned int jaccarddist = intersection_word(first->data[i], second->data[i]);

		if (jaccarddist > 0) {
			jaccarddist /= union_word(first->data[i], second->data[i]);
		}
		if (jaccarddist >= sc->thworddist) {
			if (++cwthcount > 2) {
				/* more than half (5/2) of distances are too wide */
				return 0;
			}
		}
		composdist += jaccarddist;
		if (composdist > sc->thcomposdist) {
			return 0;
		}
	}
	return 1;
}

static int find_next_coarsecandidate(SignatureContext *sc, CoarseSignature *secondstart, CoarseSignature **first, CoarseSignature **second,	int start)
{
	/* go one coarsesignature foreword */
	if (!start) {
		if ((*second)->next) {
			*second = (*second)->next;
		}
		else if ((*first)->next) {
			*second = secondstart;
			*first = (*first)->next;
		}
		else {
			return 0;
		}
	}

	while (1) {
		if (get_jaccarddist(sc, *first, *second))
			return 1;

		/* next signature */
		if ((*second)->next) {
			*second = (*second)->next;
		}
		else if ((*first)->next) {
			*second = secondstart;
			*first = (*first)->next;
		}
		else {
			return 0;
		}
	}
}

/* hough transformation */
size_t houghTransform(struct pairs *pairs, hspace_elem hspace[][2 * HOUGH_MAX_OFFSET + 1]) {
	int framerate = 0, offset = 0;
	size_t score = 0, hmax = 0;
	float m = 0;
	// This entire codeblock is a mess and it deserves a damnatio
	// memoriae. This should be a basic hough transform but for some reason
	// it's an intricate 4 level nested loop with non sequential memory
	// access and with wrong math calculations (the code tries
	// to do a linear regression by using an int division and if you
	// try to fix it it doesn't work at all). I'll try to fix it properly
	// when i'll have a clearer mind and more time.

	for (unsigned int i = 0; i < COARSE_SIZE; i++) {
		struct pairs pairI = pairs[i];

		for (unsigned int k = i + 1; k < COARSE_SIZE; k++) {
			struct pairs pairK = pairs[k];

			for (unsigned int j = 0; j < pairI.size; j++) {

				for (unsigned int l = 0; l < pairK.size; l++) {
					if (pairI.b[j] != pairK.b[l]) {

						// linear regression
						// good value between 0.0 - 2.0
						m = ((float)(pairK.b_pos[l] - pairI.b_pos[j])) / \
							(k - i);
						// round up to 0 - 60
						framerate = nearbyint(m * 30 + 0.5);
						if (framerate > 0 && framerate <= MAX_FRAMERATE) {
							// only second part has to be rounded up
							offset = pairI.b_pos[j] - nearbyint(m*i + 0.5);

							if (offset > -HOUGH_MAX_OFFSET && offset < HOUGH_MAX_OFFSET) {
								hspace_elem *hElem = &hspace[framerate - 1][offset + HOUGH_MAX_OFFSET];
								if (pairI.dist < pairK.dist) {
									if (pairI.dist < hElem->dist) {
										hElem->dist = pairI.dist;
										hElem->a = pairI.a;
										hElem->b = pairI.b[j];
									}
								}
								else {
									if (pairK.dist < hElem->dist) {
										hElem->dist = pairK.dist;
										hElem->a = pairK.a;
										hElem->b = pairK.b[l];
									}
								}

								score = ++hElem->score;
								if (score > hmax)
									hmax = score;
							}
						}
					}
				}
			}
		}
	}
	return hmax;
}

static MatchingInfo* get_matching_parameters(SignatureContext *sc, FineSignature *first, FineSignature *second)
{
	// This function must be optimized
	FineSignature *f;
	//size_t hmax = 0;
	//int l1dist;

	size_t hmax = 0;
	unsigned int l1dist = 999999999;

	MatchingInfo *cands = NULL, *c = NULL;
	struct pairs pairs[COARSE_SIZE];

	for (int i = 0; i < COARSE_SIZE; i++) {
		pairs[i].size = 0; pairs[i].dist = 999999999; pairs[i].a = NULL;
		for (int j = 0; j < COARSE_SIZE; j++) {
			pairs[i].b_pos[j] = 0;
			pairs[i].b[j] =  NULL;
		}
	}

	// initialize houghspace
	// Removed malloc/calloc to avoid useless memory allocation
	// when plenty of space is avaiable in the stack
	// designated initializer notation, removes initialization loops
	hspace_elem hspace[MAX_FRAMERATE][2 * HOUGH_MAX_OFFSET + 1];
	for (int i = 0; i < MAX_FRAMERATE; i++)	{
		for (int j = 0; j < 2 * HOUGH_MAX_OFFSET; j++) {
			hspace[i][j].score = 0;
			hspace[i][j].dist = 999999999;
			hspace[i][j].a = NULL;
			hspace[i][j].b = NULL;
		}
	}

	// l1 distances
	f = first;
	for (unsigned int i = 0; i < COARSE_SIZE && f->next; i++, f = f->next) {
		pairs[i].a = f;
		FineSignature *s = second;
		for (unsigned int j = 0; j < COARSE_SIZE && s->next; \
			j++, s = s->next) {
			// l1 distance of finesignature
			l1dist = (unsigned int)get_l1dist(sc, f->framesig, s->framesig);
			if (l1dist < (unsigned int)sc->thl1) {
				if (l1dist < pairs[i].dist) {
					pairs[i].size = 1;
					pairs[i].dist = l1dist;
					pairs[i].b_pos[0] = j;
					pairs[i].b[0] = s;
				}
				else if (l1dist == pairs[i].dist) {
					pairs[i].b[pairs[i].size] = s;
					pairs[i].b_pos[pairs[i].size] = j;
					pairs[i].size++;
				}
			}
		}
	}

	hmax = houghTransform(pairs, hspace);

	if (hmax > 0) {
		hmax = floor(0.7*hmax);
		for (unsigned int i = 0; i < MAX_FRAMERATE; i++) {
			for (unsigned int j = 0; j < HOUGH_MAX_OFFSET; j++) {
				if (hmax < hspace[i][j].score) {
					if (c == NULL) {
						c = (MatchingInfo*)calloc(1, sizeof(MatchingInfo));						
						cands = c;
					}
					else {
						c->next = (MatchingInfo*)calloc(1, sizeof(MatchingInfo));						
						c = c->next;
					}
					c->framerateratio = (i + 1.0) / 30;
					c->score = hspace[i][j].score;
					c->offset = j - 90;
					c->first = hspace[i][j].a;
					c->second = hspace[i][j].b;
					c->next = NULL;

					// not used
					c->meandist = 0;
					c->matchframes = 0;
					c->whole = 0;
				}
			}
		}
	}
	return cands;
}

static int iterate_frame(double frr, FineSignature **a,	FineSignature **b, int fcount, int *bcount,	int dir)
{
	int step;

	/* between 1 and 2, because frr is between 1 and 2 */
	step = ((int) 0.5 + fcount * frr) /* current frame */
		- ((int) 0.5 + (fcount - 1) * frr);/* last frame */

	if (dir == DIR_NEXT) {
		if (frr >= 1.0) {
			if ((*a)->next) {
				*a = (*a)->next;
			}
			else {
				return DIR_NEXT_END;
			}

			if (step == 1) {
				if ((*b)->next) {
					*b = (*b)->next;
					(*bcount)++;
				}
				else {
					return DIR_NEXT_END;
				}
			}
			else {
				if ((*b)->next && (*b)->next->next) {
					*b = (*b)->next->next;
					(*bcount)++;
				}
				else {
					return DIR_NEXT_END;
				}
			}
		}
		else {
			if ((*b)->next) {
				*b = (*b)->next;
				(*bcount)++;
			}
			else {
				return DIR_NEXT_END;
			}

			if (step == 1) {
				if ((*a)->next) {
					*a = (*a)->next;
				}
				else {
					return DIR_NEXT_END;
				}
			}
			else {
				if ((*a)->next && (*a)->next->next) {
					*a = (*a)->next->next;
				}
				else {
					return DIR_NEXT_END;
				}
			}
		}
		return DIR_NEXT;
	}
	else {
		if (frr >= 1.0) {
			if ((*a)->prev) {
				*a = (*a)->prev;
			}
			else {
				return DIR_PREV_END;
			}

			if (step == 1) {
				if ((*b)->prev) {
					*b = (*b)->prev;
					(*bcount)++;
				}
				else {
					return DIR_PREV_END;
				}
			}
			else {
				if ((*b)->prev && (*b)->prev->prev) {
					*b = (*b)->prev->prev;
					(*bcount)++;
				}
				else {
					return DIR_PREV_END;
				}
			}
		}
		else {
			if ((*b)->prev) {
				*b = (*b)->prev;
				(*bcount)++;
			}
			else {
				return DIR_PREV_END;
			}

			if (step == 1) {
				if ((*a)->prev) {
					*a = (*a)->prev;
				}
				else {
					return DIR_PREV_END;
				}
			}
			else {
				if ((*a)->prev && (*a)->prev->prev) {
					*a = (*a)->prev->prev;
				}
				else {
					return DIR_PREV_END;
				}
			}
		}
		return DIR_PREV;
	}
}


static MatchingInfo evaluate_parameters(SignatureContext *sc, MatchingInfo *infos, MatchingInfo bestmatch)
{
	int dist, distsum = 0, bcount = 1, dir = DIR_NEXT;
	int fcount = 0, goodfcount = 0, gooda = 0, goodb = 0;
	double meandist, minmeandist = bestmatch.meandist;
	int tolerancecount = 0;
	FineSignature *aprev = NULL, *bprev = NULL;
	int status = STATUS_NULL;

	for (; infos != NULL; infos = infos->next) {
		FineSignature *a = infos->first;
		FineSignature *b = infos->second;
		while (1) {
			dist = get_l1dist(sc, a->framesig, b->framesig);

			if (dist > sc->thl1) {
				if (a->confidence >= 1 || b->confidence >= 1) {
					/* bad frame (because high different information) */
					tolerancecount++;
				}

				if (tolerancecount > 2) {
					a = aprev;
					b = bprev;
					if (dir == DIR_NEXT) {
						/* turn around */
						a = infos->first;
						b = infos->second;
						dir = DIR_PREV;
					}
					else {
						break;
					}
				}
			}
			else {
				/* good frame */
				distsum += dist;
				goodfcount++;
				tolerancecount = 0;

				aprev = a;
				bprev = b;

				if (a->confidence < 1) gooda++;
				if (b->confidence < 1) goodb++;
			}

			fcount++;

			dir = iterate_frame(infos->framerateratio, &a, &b, fcount, &bcount, dir);
			if (dir == DIR_NEXT_END) {
				status = STATUS_END_REACHED;
				a = infos->first;
				b = infos->second;
				dir = iterate_frame(infos->framerateratio, &a, &b, fcount, &bcount, DIR_PREV);
			}

			if (dir == DIR_PREV_END) {
				status |= STATUS_BEGIN_REACHED;
				break;
			}

			if (sc->thdi != 0 && bcount >= sc->thdi) {
				break; /* enough frames found */
			}
		}

		if (bcount < sc->thdi)
			continue; /* matching sequence is too short */
		if ((double)goodfcount / (double)fcount < sc->thit)
			continue;
		if ((double)goodfcount*0.5 < FFMAX(gooda, goodb))
			continue;

		meandist = (double)goodfcount / (double)distsum;

		if (meandist < minmeandist ||
			status == (STATUS_END_REACHED | STATUS_BEGIN_REACHED) ||
			(sc->mode == MODE_FAST)) {
			minmeandist = meandist;
			/* bestcandidate in this iteration */
			bestmatch.meandist = meandist;
			bestmatch.matchframes = bcount;
			bestmatch.framerateratio = infos->framerateratio;
			bestmatch.score = infos->score;
			bestmatch.offset = infos->offset;
			bestmatch.first = infos->first;
			bestmatch.second = infos->second;
			bestmatch.whole = 0; /* will be set to true later */
			bestmatch.next = NULL;
		}

		/* whole sequence is automatically best match */
		if (status == (STATUS_END_REACHED | STATUS_BEGIN_REACHED)) {
			bestmatch.whole = 1;
			break;
		}

		/* first matching sequence is enough, finding the best one is not necessary */
		if (sc->mode == MODE_FAST) {
			break;
		}
	}
	return bestmatch;
}

static void sll_free(MatchingInfo *sll)
{
	while (sll) {
		void *tmp = sll;
		sll = sll->next;
		free(tmp);
	}
}

MatchingInfo lookup_signatures(SignatureContext *sc, StreamContext *first, StreamContext *second)
{
	CoarseSignature *cs = NULL, *cs2 = NULL;
	MatchingInfo *infos = NULL;
	MatchingInfo bestmatch = { 0 };

	cs = first->coarsesiglist;
	cs2 = second->coarsesiglist;

	/* score of bestmatch is 0, if no match is found */
	bestmatch.score = 0;
	bestmatch.meandist = 99999;
	bestmatch.whole = 0;

	fill_l1distlut(sc->l1distlut);

	/* stage 1: coarsesignature matching */
	if (!find_next_coarsecandidate(sc, second->coarsesiglist, &cs, &cs2, 1))
		return bestmatch; /* no candidate found */

	do {
		// stage 2: l1-distance and hough-transform
		infos = get_matching_parameters(sc, cs->first, cs2->first);
		
		// stage 3: evaluation
		if (infos) {
			bestmatch = evaluate_parameters(sc, infos, bestmatch);		
			sll_free(infos);
		}
	} while (find_next_coarsecandidate(sc, second->coarsesiglist, &cs, &cs2, 0) && !bestmatch.whole);

	return bestmatch;
}
