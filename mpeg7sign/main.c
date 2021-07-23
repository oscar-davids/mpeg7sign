// mpeg7sign.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include "libavcodec/avcodec.h"
#include "signature.h"
#include "get_bits.h"

#define INPUTS_COUNT  2

typedef struct BoundedCoarseSignature {
	// StartFrameOfSegment and EndFrameOfSegment
	uint32_t firstIndex, lastIndex;
	// StartMediaTimeOfSegment and EndMediaTimeOfSegment
	uint64_t firstPts, lastPts;
	CoarseSignature *cSign;
} BoundedCoarseSignature;

int load_Signaturebin(char* filename, StreamContext *sc);
void release_StreamContext(StreamContext *sc);

void release_StreamContext(StreamContext *sc)
{
	free(sc->coarsesiglist);
	free(sc->finesiglist);
}


int main(int argc, char **argv)
{
	StreamContext scontexts[INPUTS_COUNT] = { 0 };
	MatchingInfo result = { 0 };
	SignatureContext sigContext = {
		.class = NULL,
		.mode = MODE_FULL,
		.nb_inputs = INPUTS_COUNT,
		.filename = "",
		.thworddist = 9000,
		.thcomposdist = 60000,
		.thl1 = 116,
		.thdi = 0,
		.thit = 0.5,
		.streamcontexts = scontexts
	};

	if (load_Signaturebin(argv[1], &scontexts[0]) < 0 || load_Signaturebin(argv[2], &scontexts[1]) < 0) {
		printf("Could not load signature binary file!\n");
	}

	result = lookup_signatures(&sigContext, &scontexts[0], &scontexts[1]);

	if (result.score != 0) {
		if (result.whole) {
			printf("whole matching\n");
		}
		else {
			printf("matching frames %d\n", result.matchframes);
		}
	}
	else {
		printf("no matching\n");
	}	

	release_StreamContext(&scontexts[0]);
	release_StreamContext(&scontexts[1]);
	return 0;
}

int fineSignatureCmp(const void* p1, const void* p2) {
	FineSignature *a = (FineSignature*)p1;
	FineSignature *b = (FineSignature*)p2;
	if (a->pts == b->pts)
		return 0;
	else if (a->pts < b->pts)
		return -1;
	else
		return 1;
}

unsigned int getFileSize(const char *filename) 
{
	int fileLength = 0;
	FILE *f = NULL;
	f = fopen(filename, "rb");	
	fseek(f, 0, SEEK_END);
	fileLength = ftell(f);
	fclose(f);
	return fileLength;
};

int load_Signaturebin(char* filename, StreamContext *sc)
{
	int ret = 0;

	FILE *f = NULL;
	unsigned int rResult = 0, fileLength = 0, paddedLength = 0, \
		numOfSegments = 0;
	uint8_t *buffer = NULL;
	GetBitContext bitContext = { 0 };
	//check input parameters
	if (strlen(filename) <= 0 || sc == NULL) return -1;
	f = fopen(filename, "rb");
	if (f == NULL) return -1;

	fileLength = getFileSize(filename);

	// Cast to float is necessary to avoid int division
	paddedLength = ceil(fileLength / (float)AV_INPUT_BUFFER_PADDING_SIZE)*AV_INPUT_BUFFER_PADDING_SIZE + AV_INPUT_BUFFER_PADDING_SIZE;
	buffer = (uint8_t*)calloc(paddedLength, sizeof(uint8_t));	

	// Read entire file into memory
	rResult = fread(buffer, sizeof(uint8_t), fileLength, f);
	//Assert(rResult == fileLength);
	// Remove FILE pointer from memory once we're done
	fclose(f);
	f = NULL;

	// BE CAREFUL, THE LENGTH IS SPECIFIED IN BITS NOT BYTES
	init_get_bits(&bitContext, buffer, 8 * fileLength);
	// libavcodec

	// Skip the following data:
	// - NumOfSpatial Regions: (32 bits) only 1 supported
	// - SpatialLocationFlag: (1 bit) always the whole image
	// - PixelX_1: (16 bits) always 0
	// - PixelY_1: (16 bits) always 0
	skip_bits(&bitContext, 32 + 1 + 16 * 2);


	// width - 1, and height - 1
	// PixelX_2: (16 bits) is width - 1
	// PixelY_2: (16 bits) is height - 1
	sc->w = get_bits(&bitContext, 16);
	sc->h = get_bits(&bitContext, 16);
	++sc->w;
	++sc->h;

	// StartFrameOfSpatialRegion, always 0
	skip_bits(&bitContext, 32);

	// NumOfFrames
	// it's the number of fine signatures
	sc->lastindex = get_bits_long(&bitContext, 32);

	// sc->time_base.den / sc->time_base.num
	// hoping num is 1, other values are vague
	// den/num might be greater than 16 bit, so cutting it
	//put_bits(&buf, 16, 0xFFFF & (sc->time_base.den / sc->time_base.num));
	// MediaTimeUnit

	sc->time_base.den = get_bits(&bitContext, 16);
	sc->time_base.num = 1;

	// Skip the following data
	// - MediaTimeFlagOfSpatialRegion: (1 bit) always 1
	// - StartMediaTimeOfSpatialRegion: (32 bits) always 0
	skip_bits(&bitContext, 1 + 32);

	// EndMediaTimeOfSpatialRegion
	uint64_t lastCoarsePts = get_bits_long(&bitContext, 32);

	// Coarse signatures
	// numOfSegments = number of coarse signatures
	numOfSegments = get_bits_long(&bitContext, 32);

	// Reading from binary signature return a wrong number of segments
	// numOfSegments = (sc->lastindex + 44)/45;
	//skip_bits(&bitContext, 32);

	sc->coarsesiglist = (CoarseSignature*)calloc(numOfSegments,	sizeof(CoarseSignature));

	BoundedCoarseSignature *bCoarseList = (BoundedCoarseSignature*)calloc(numOfSegments, sizeof(BoundedCoarseSignature));


	// CoarseSignature loading	

	for (unsigned int i = 0; i < numOfSegments; ++i) {
		BoundedCoarseSignature *bCs = &bCoarseList[i];
		bCs->cSign = &sc->coarsesiglist[i];

		if (i < numOfSegments - 1)
			bCs->cSign->next = &sc->coarsesiglist[i + 1];

		// each coarse signature is a VSVideoSegment

		// StartFrameOfSegment
		bCs->firstIndex = get_bits_long(&bitContext, 32);
		// EndFrameOfSegment
		bCs->lastIndex = get_bits_long(&bitContext, 32);

		// MediaTimeFlagOfSegment 1 bit, always 1
		skip_bits(&bitContext, 1);

		// Fine signature pts
		// StartMediaTimeOfSegment 32 bits
		bCs->firstPts = get_bits_long(&bitContext, 32);
		// EndMediaTimeOfSegment 32 bits
		bCs->lastPts = get_bits_long(&bitContext, 32);


		// Bag of words
		for (unsigned int i = 0; i < 5; ++i) {

			// read 243 bits ( = 7 * 32 + 19 = 8 * 28 + 19) into buffer
			for (unsigned int j = 0; j < 30; ++j) {
				// 30*8 bits = 30 bytes
				bCs->cSign->data[i][j] = get_bits(&bitContext, 8);
			}
			bCs->cSign->data[i][30] = get_bits(&bitContext, 3) << 5;
		}
	}
	sc->coarseend = &sc->coarsesiglist[numOfSegments - 1];

	// Finesignatures
	// CompressionFlag, only 0 supported
	skip_bits(&bitContext, 1);


	sc->finesiglist = (FineSignature*)calloc(sc->lastindex,	sizeof(FineSignature));	

	// Load fine signatures from file	

	for (unsigned int i = 0; i < sc->lastindex; ++i) {
		FineSignature *fs = &sc->finesiglist[i];

		// MediaTimeFlagOfFrame always 1
		skip_bits(&bitContext, 1);

		// MediaTimeOfFrame (PTS)
		fs->pts = get_bits_long(&bitContext, 32);

		// FrameConfidence
		fs->confidence = get_bits(&bitContext, 8);

		// words
		for (unsigned int l = 0; l < 5; l++) {
			fs->words[l] = get_bits(&bitContext, 8);
		}

		// Crashes for some signature, it's a memory adding problems
		// framesignature
		for (unsigned int l = 0; l < SIGELEM_SIZE / 5; l++) {
			fs->framesig[l] = get_bits(&bitContext, 8);
		}
	};

	// Sort by frame time (pts)
	qsort(sc->finesiglist, sc->lastindex, sizeof(FineSignature), \
		fineSignatureCmp);

	// Creating FineSignature linked list
	for (unsigned int i = 0; i < sc->lastindex; ++i) {
		FineSignature *fs = &sc->finesiglist[i];
		// Building fine signature list
		// First element prev should be NULL
		// Last element next should be NULL
		if (i == 0) {
			fs->next = &fs[1];
			fs->prev = NULL;
		}
		else if (i == sc->lastindex - 1) {
			fs->next = NULL;
			fs->prev = &fs[-1];
		}
		else {
			fs->next = &fs[1];
			fs->prev = &fs[-1];
		}
	}

	// Fine signature ranges DO overlap
	// Assign FineSignatures to CoarseSignature s
	for (unsigned int i = 0; i < numOfSegments; ++i) {
		BoundedCoarseSignature *bCs = &bCoarseList[i];

		// O = n^2 probably it can be done faster
		for (unsigned int j = 0;  j < sc->lastindex  &&\
			sc->finesiglist[j].pts <= bCs->lastPts; ++j) {
			FineSignature *fs = &sc->finesiglist[j];

			if (fs->pts >= bCs->firstPts) {
				// Check if the fragment's pts is inside coarse signature
				// bounds. Upper bound is checked in for loop
				if (!bCs->cSign->first) {
					bCs->cSign->first = fs;
				}

				if (bCs->cSign->last) {
					if (bCs->cSign->last->pts <= fs->pts)
						bCs->cSign->last = fs;
				}
				else {
					bCs->cSign->last = fs;
				}
			}

		}
		bCs->cSign->first->index = bCs->firstIndex;
		bCs->cSign->last->index = bCs->lastIndex;
	}
	
	free(bCoarseList);
	free(buffer);


	return ret;
}
