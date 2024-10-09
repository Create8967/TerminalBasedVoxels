#include <fcntl.h>
#include <unistd.h>
#include <termios.h>
#include <poll.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define max(x,m) ((x>m)?x:m)

// Voxel settings
#define MAXDIST 80.0f
#define VOXELSIZE 1.0f
#define CHUNKSIZE 16
#define AIR 16

// Use less division
float divChunksize = 1.0f/(float)CHUNKSIZE;
const char shades16[16] = {' ','.',',',':','-','=','+','%','|','!','?','L','Z','&','#','@'};

int Width = 170;
int Height = 48;
int MAXSTEP = 31;
int CHUNKRAD = 1;

int world_seed = 2853;

unsigned char *pixels;

// Set cursor position
void setpos(int x, int y) { printf("\033[%d;%dH",y,x); }
void drawChar(unsigned char col, unsigned char shade) {
  if(col > 15 || shade > 15) {
  	printf("~");
  	return;
  }
  int c = (col < 8) ? 30+col : 90+(col-8);
  printf("\033[0;%dm%c",c,shades16[shade]);
}
void draw(int x, int y, unsigned char col, unsigned char shade) {
	if(x<0||x>=Width) return;
	if(y<0||y>=Height) return;
	if(shade > 15) shade = 15;
	if(shade < 0) shade = 0;
	col = col%16;

	// Color stored in top 4 bits, shade in the lower 4
	pixels[y*Width+x] = (col&15)<<4|(shade&15);
}
void clear() {
	for(int i=0;i<Width*Height;i++) {
		pixels[i] = 0;
	}
}

// Update all text pixels
void update() {
	for(int h=0;h<Height;h++) {
		setpos(0,h);
		for(int w=0;w<Width;w++) {
			unsigned char c = pixels[h*Width+w];
			drawChar((c&(15<<4))>>4, c&15);
		}
	}
	printf("\033[0m");
}

// Vector functions

typedef struct {
	float x,y,z;
} vec3;

typedef struct {
	float x,y;
} vec2;

vec3 vec3_add(vec3 a, vec3 b) {
	return (vec3){a.x+b.x,a.y+b.y,a.z+b.z};
}
vec3 vec3_sub(vec3 a, vec3 b) {
	return (vec3){a.x-b.x,a.y-b.y,a.z-b.z};
}
vec3 vec3_mul(vec3 a, vec3 b) {
	return (vec3){a.x*b.x,a.y*b.y,a.z*b.z};
}
vec3 vec3_div(vec3 a, vec3 b) {
	return (vec3){a.x/b.x,a.y/b.y,a.z/b.z};
}
vec3 vec3_mulf(vec3 a, float b) {
	return (vec3){a.x*b,a.y*b,a.z*b};
}
vec3 vec3_divf(vec3 a, float b) {
	return (vec3){a.x/b,a.y/b,a.z/b};
}
vec2 vec2_sub(vec2 a, vec2 b) {
	return (vec2){a.x-b.x,a.y-b.y};
}
float length(vec3 a) {
	return sqrtf(a.x*a.x+a.y*a.y+a.z*a.z);
}
float dot(vec3 a, vec3 b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
float dot2(vec2 a, vec2 b) {
	return a.x+b.x+a.y+b.y;
}
vec3 normalize(vec3 a) {
	float l = length(a);
	return vec3_divf(a,l);
}
vec2 Rot(float x, float y, float a) {
	float s = sin(a);
	float c = cos(a);

	return (vec2){x * c + y * s, x * -s + y * c};
}

//  ****************************


// Perlin Noise functions

float fract(float x) {
	return x - floorf(x);
}

// 2d random
vec2 random2(vec2 v, int seed) {
	float r = dot2(v,(vec2){93.1057f,37.9163f});
	float addr = *(float*)(&seed);
	float ran = fract(sinf(r+addr)*438243.6471f);
	return (vec2){cosf(ran*3.14159f*2.0f),sinf(ran*3.14159f*2.0f)};
}

// 3d random
vec3 random3(vec3 v, int seed) {
	float r = dot(v,(vec3){93.1057f,37.9163f, 21.7352f});
	float addr = *(float*)(&seed);
	float ran = fract(sinf(r+addr)*438243.6471f);
	float ran2 = fract(cosf(r+addr)*156243.6471f);
	return (vec3){cosf(ran*3.14159f*2.0f),sinf(ran*3.14159f*2.0f),cosf(ran2*3.14159f*2.0f)};
}

float lerp(float a, float b, float t) {
	return a+(b-a)*t;
}
float slerp(float a, float b, float t) {
	float x = t*t*(-2.0f*t + 3);
	return lerp(a,b,x);
}

// 2d perlin noise
float perlin_noise(float x, float y, int seed) {
	vec2 g = (vec2){x,y};

	vec2 g0 = (vec2){floorf(x),floorf(y)};
	vec2 g1 = (vec2){floorf(x)+1.0f,floorf(y)};
	vec2 g2 = (vec2){floorf(x),floorf(y)+1.0f};
	vec2 g3 = (vec2){floorf(x)+1.0f,floorf(y)+1.0f};

	vec2 r0 = random2(g0, seed);
	vec2 r1 = random2(g1, seed);
	vec2 r2 = random2(g2, seed);
	vec2 r3 = random2(g3, seed);

	float d0 = dot2(r0, vec2_sub(g,g0));
	float d1 = dot2(r1, vec2_sub(g,g1));
	float d2 = dot2(r2, vec2_sub(g,g2));
	float d3 = dot2(r3, vec2_sub(g,g3));

	float h1 = slerp(d0,d1,g.x-g0.x);
	float h2 = slerp(d2,d3,g.x-g0.x);
	float v = slerp(h1,h2,g.y-g0.y);

	return v;
}

// 3d perlin noise
float perlin_noise3(float x, float y, float z, int seed) {
	vec3 g = (vec3){x,y,z};

	vec3 g0 = (vec3){floorf(x),floorf(y),floorf(z)};
	vec3 g1 = (vec3){floorf(x)+1.0f,floorf(y),floorf(z)};
	vec3 g2 = (vec3){floorf(x),floorf(y)+1.0f,floorf(z)};
	vec3 g3 = (vec3){floorf(x)+1.0f,floorf(y)+1.0f,floorf(z)};

	vec3 g4 = (vec3){floorf(x),floorf(y),floorf(z)+1.0f};
	vec3 g5 = (vec3){floorf(x)+1.0f,floorf(y),floorf(z)+1.0f};
	vec3 g6 = (vec3){floorf(x),floorf(y)+1.0f,floorf(z)+1.0f};
	vec3 g7 = (vec3){floorf(x)+1.0f,floorf(y)+1.0f,floorf(z)+1.0f};

	vec3 r0 = random3(g0, seed);
	vec3 r1 = random3(g1, seed);
	vec3 r2 = random3(g2, seed);
	vec3 r3 = random3(g3, seed);

	vec3 r4 = random3(g4, seed);
	vec3 r5 = random3(g5, seed);
	vec3 r6 = random3(g6, seed);
	vec3 r7 = random3(g7, seed);

	float d0 = dot(r0, vec3_sub(g,g0));
	float d1 = dot(r1, vec3_sub(g,g1));
	float d2 = dot(r2, vec3_sub(g,g2));
	float d3 = dot(r3, vec3_sub(g,g3));

	float d4 = dot(r4, vec3_sub(g,g4));
	float d5 = dot(r5, vec3_sub(g,g5));
	float d6 = dot(r6, vec3_sub(g,g6));
	float d7 = dot(r7, vec3_sub(g,g7));

	float h1 = slerp(d0,d1,g.x-g0.x);
	float h2 = slerp(d2,d3,g.x-g0.x);
	float v1 = slerp(h1,h2,g.y-g0.y);

	float h3 = slerp(d4,d5,g.x-g0.x);
	float h4 = slerp(d6,d7,g.x-g0.x);
	float v2 = slerp(h3,h4,g.y-g0.y);

	float c = slerp(v1,v2,g.z-g0.z);

	return c;
}

// Octave noise
float oct_noise(float x, float y, int octaves, float scale, float perm, int seed) {
	float fscale = scale;
	float fperm = perm;
	float max_scale = 1.0f;
	float noise = perlin_noise(x,y,seed);
	for(int i=1;i<octaves;i++) {
		noise += perlin_noise(x*fscale,y*fscale,seed) * fperm;

		max_scale += fperm;
		fscale *= scale;
		fperm *= perm;
	}
	return noise / max_scale;
}

// 3d octave noise
float oct_noise3(float x, float y, float z, int octaves, float scale, float perm, int seed) {
	float fscale = scale;
	float fperm = perm;
	float max_scale = 1.0f;
	float noise = perlin_noise3(x,y,z,seed);
	for(int i=1;i<octaves;i++) {
		noise += perlin_noise3(x*fscale,y*fscale,z*fscale,seed) * fperm;

		max_scale += fperm;
		fscale *= scale;
		fperm *= perm;
	}
	return noise / max_scale;
}

// Voxel Stuff

typedef struct {
	char *voxel;
	int x,y,z;
} Chunk;

// Chunk creator
int newChunk(int cx, int cy, int cz, Chunk *newC) {
	newC->x=cx; newC->y=cy; newC->z=cz;
	newC->voxel = (char*)calloc(CHUNKSIZE*CHUNKSIZE*CHUNKSIZE,sizeof(char));
	if(newC->voxel == NULL) return -1;
	return 0;
}
int chunkPos(int x, int y, int z) {
	return x + y*CHUNKSIZE + z*CHUNKSIZE*CHUNKSIZE;
}
int getBlock(float x, float y, float z, Chunk *wchunks, int wchunksize, Chunk *pchunks, int pchunksize) {
	int vx = floorf(x/VOXELSIZE);
	int vy = floorf(y/VOXELSIZE);
	int vz = floorf(z/VOXELSIZE);
	int cx = floorf((float)vx*divChunksize);
	int cy = floorf((float)vy*divChunksize);
	int cz = floorf((float)vz*divChunksize);
	int wwchunk = -1;
	int ppchunk = -1;
	for(int e=0;e<wchunksize;e++) {
		if((cx == wchunks[e].x && cy == wchunks[e].y) && cz == wchunks[e].z) {
			wwchunk = e;
			break;
		}
	}
	for(int p=0;p<pchunksize;p++) {
		if((cx == pchunks[p].x && cy == pchunks[p].y) && cz == pchunks[p].z) {
			ppchunk = p;
			break;
		}
	}
	if(wwchunk == -1) return 0; // L bozo
	else {
		if(ppchunk != -1) {
			int tmpvoxel = pchunks[ppchunk].voxel[chunkPos((vx%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vy%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vz%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE)];
			return (tmpvoxel != 0) ? tmpvoxel : wchunks[wwchunk].voxel[chunkPos((vx%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vy%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vz%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE)];
		}
		else {
			return wchunks[wwchunk].voxel[chunkPos((vx%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vy%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vz%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE)];
		}
	}
}

// Test chunk function
void loadSphere(int r, Chunk *c) {
	for(int z= -r; z<=r;z++) {
		for(int y= -r; y<=r;y++) {
			for(int x= -r; x<=r;x++) {
				if(x*x+y*y+z*z < r*r) {
					c->voxel[chunkPos(x+CHUNKSIZE/2,y+CHUNKSIZE/2,z+CHUNKSIZE/2)] = 1;
				}
			}
		}
	}
}

// Voxel raycasting using "A Fast Voxel Traversal Algorithm for Ray Tracing" by John Amanatides and Andrew Woo
float voxelTrace(vec3 ro, vec3 rd, Chunk *worldChunks, int wChunkSize, Chunk *playerChunks, int pChunkSize, vec3 *normal, int *vox, int *vx, int *vy, int *vz) {

	// Current voxel position
	int voxX = floorf(ro.x/VOXELSIZE);
	int voxY = floorf(ro.y/VOXELSIZE);
	int voxZ = floorf(ro.z/VOXELSIZE);

	// Voxel ray direction
	int stepX = (rd.x < 0.0f)? -1:1;
	int stepY = (rd.y < 0.0f)? -1:1;
	int stepZ = (rd.z < 0.0f)? -1:1;

	int wchunk = -1; // current world chunk
	int pchunk = -1; // current player chunk

	// Chunk pos
	int pcx = voxX/CHUNKSIZE-2;
	int pcy = voxY/CHUNKSIZE-2;
	int pcz = voxZ/CHUNKSIZE-2;

	// Minimum distance to next voxel
	vec3 tMax = (vec3){(stepX == 1)? (float)(voxX+1)*VOXELSIZE - ro.x : (float)voxX*VOXELSIZE - ro.x,
											 (stepY == 1)? (float)(voxY+1)*VOXELSIZE - ro.y : (float)voxY*VOXELSIZE - ro.y,
											 (stepZ == 1)? (float)(voxZ+1)*VOXELSIZE - ro.z : (float)voxZ*VOXELSIZE - ro.z};
	tMax = vec3_div(tMax,rd);

	// Distance for one voxel unit
	vec3 tDelta = (vec3){VOXELSIZE/fabs(rd.x), VOXELSIZE/fabs(rd.y), VOXELSIZE/fabs(rd.z)};

	int voxel = AIR;
	float dist = 0.0f;
	vec3 n = {0.0f, 0.0f, 0.0f};

	// Step until you land on something that isn't AIR
	for(int i=0; i<MAXSTEP;i++) {
		if(tMax.x < tMax.y) {
			if(tMax.x < tMax.z) {
				dist = tMax.x;
				tMax.x += tDelta.x;
				voxX += stepX;
				n = (vec3){-stepX, 0, 0};
			}
			else {
				dist = tMax.z;
				tMax.z += tDelta.z;
				voxZ += stepZ;
				n = (vec3){0,0, -stepZ};
			}
		}
		else {
			if(tMax.y < tMax.z) {
				dist = tMax.y;
				tMax.y += tDelta.y;
				voxY += stepY;
				n = (vec3){0, -stepY, 0};
			}
			else {
				dist = tMax.z;
				tMax.z += tDelta.z;
				voxZ += stepZ;
				n = (vec3){0,0, -stepZ};
			}
		}

		// Check current voxel
		int cx = floorf((float)voxX*divChunksize);
		int cy = floorf((float)voxY*divChunksize);
		int cz = floorf((float)voxZ*divChunksize);

		// If you've changed chunks
		if((cx!=pcx || cy != pcy)||cz != pcz) {
			wchunk = -1;
			pchunk = -1;
			for(int j = 0; j<wChunkSize;j++) {
				if((cx==worldChunks[j].x && cy==worldChunks[j].y) && cz==worldChunks[j].z) {
					wchunk = j;
					break;
				}
			}
			for(int k = 0; k<pChunkSize;k++) {
				if((cx==playerChunks[k].x && cy==playerChunks[k].y) && cz==playerChunks[k].z) {
					pchunk = k;
					break;
				}
			}
			pcx = cx;
			pcy = cy;
			pcz = cz;
		}
		if(wchunk == -1) { // Shouldn't happen but if it does dont render allat
			voxel = AIR;
			break;
		}
		else {
			if(pchunk != -1) { // Player chunks come before world chunks
				int tempV = playerChunks[pchunk].voxel[chunkPos((voxX%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE, (voxY%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE, (voxZ%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE)];
				voxel = (tempV != 0)? tempV : worldChunks[wchunk].voxel[chunkPos((voxX%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE, (voxY%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE, (voxZ%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE)];
			}
			else { // World voxel
				voxel = worldChunks[wchunk].voxel[chunkPos((voxX%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE, (voxY%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE, (voxZ%CHUNKSIZE + CHUNKSIZE)%CHUNKSIZE)];
			}
		}
		if(voxel != AIR) { // FINALLY hit something
			break;
		}
	}
	*vox = voxel;
	*normal = n;
	*vx = voxX;
	*vy = voxY;
	*vz = voxZ;
	if(voxel == AIR) dist = MAXDIST;
	return dist;
}


float getLight(vec3 n, vec3 l) {
	float lig = dot(n,l) * 0.5f + 0.5f;
	return lig;
}

// **** World Types ****

char flatWorld(int x, int y, int z) {
	if(y<=0) return 1;
	return AIR;
}
char trigWorld(int x, int y, int z) {
	int h = 3.0f*sin((float)x*0.05f)+2.0f*cos((float)z*0.1f) + 3;
	if(h<3) h=3;
	if(y<=h) {
		return ((y)%16+16)%16;
	}
	return AIR;
}
char hillyWorld(int x, int y, int z) {
	float dx = (float)x * 0.05f;
	float dy = (float)y * 0.05f;
	float dz = (float)z * 0.05f;
	int height = 0;
	float caves = 0.0f;
	float clouds = 0.0f;
	if(y> -4) {
		float hills = max(oct_noise(dx,dz, 1, 2.0f, 0.5f, world_seed), -0.1f);
		float mountains = -(oct_noise(dx*0.1f,dz*0.1f, 6, 2.0f, 0.6f, world_seed+127));
		float mountain_clamp = 0.1f;
		if(mountains < mountain_clamp) mountains = 0.0f;
		else mountains = (mountains-mountain_clamp)*(1.0f/(1.0f-mountain_clamp));

		height = hills*4.0f + 1 + (mountains*mountains*mountains)*150.0f;

	}
	else {
		caves = oct_noise3(dx,dy,dz, 2, 2.0f, 0.5f, world_seed+1534);
	}
	if(y>18) {
		clouds = oct_noise3(dx,dy,dz, 2, 2.0f, 0.5f, world_seed+2351);
	}
	if(y > height+20) {
		if(clouds>0.1f) return (int)(clouds*10.0f)+6;
		return AIR;
	}
	else if(y<=height) {
		if(y==height) {
			return 2; // Grass
		}
		else if(y < height && y > height -2) {
			return 3; // Dirt
		}
		else if(y< -5) {
			if(caves > 0.1f) {
				return AIR;
			}
		}
		return 8; // Stone
	}
	return AIR;
}

// Test function for 2d height map
void renderMap(int scale, int move) {
	for(int y = 0; y<Height;y++) {
		for(int x = 0; x<Width;x++) {
			float dx = (float)(x*scale+move) * 0.05f;
			float dz = (float)(y*scale) * 0.05f;
			float mountains = -(oct_noise(dx*0.4f,dz*0.4f, 3, 2.0f, 0.6f, world_seed+127));
			float mountain_clamp = 0.4f;
			if(mountains < mountain_clamp) mountains = 0.0f;
			else mountains = (mountains-mountain_clamp)*(1.0f/(1.0f-mountain_clamp));
			int shad = mountains*15.0f;
			if(shad > 15) shad = 15;
			if(shad < 0) shad = 0;
			draw(x,y,15,shad);
		}
	}
}

// Load world chunk from world type function
Chunk loadChunkWorld(int cx, int cy, int cz, char worldType(int,int,int)) {
	Chunk c;
	int status = newChunk(cx,cy,cz,&c);
	if(status == -1) exit(0);
	for(int z=0;z<CHUNKSIZE;z++) {
		for(int y=0;y<CHUNKSIZE;y++) {
			for(int x=0;x<CHUNKSIZE;x++) {
				c.voxel[chunkPos(x,y,z)]=worldType(cx*CHUNKSIZE+x,cy*CHUNKSIZE+y,cz*CHUNKSIZE+z);
			}
		}
	}
	return c;
}

// Add chunk to chunk array
Chunk *addChunk(int cx,int cy,int cz, Chunk *c, int *chunksize, int *chunkmax) {
	int chunkRad = MAXSTEP/(VOXELSIZE*CHUNKSIZE);
	int chunkMaxima = (2*chunkRad+1)*(2*chunkRad+1)*(2*chunkRad+1);
	Chunk *nc = c;
	int newchunksize = *chunksize;
	if(newchunksize >= *chunkmax) { // gotta reallac
		*chunkmax += chunkMaxima;
		nc = (Chunk*)realloc(c, (*chunkmax)*sizeof(Chunk));
	}
	Chunk cccc;
	if(newChunk(cx,cy,cz,nc+(newchunksize++))==-1) {
		printf("YOU FAILED AT LIFE\n");
	}
	*chunksize = newchunksize;
	return nc;
}

void setBlock(char block, int vx, int vy, int vz, Chunk **c, int *chunksize, int *chunkmax) {
	// Check for chunk
	int cx = floorf((float)vx*divChunksize);
	int cy = floorf((float)vy*divChunksize);
	int cz = floorf((float)vz*divChunksize);
	int cchunk = -1;
	for(int e=0;e<(*chunksize);e++) {
		if(((*c+e)->x == cx && (*c+e)->y==cy)&& (*c+e)->z == cz) {
			cchunk = e;
			break;
		}
	}
	if(cchunk==-1) {
		*c = addChunk(cx,cy,cz,*c,chunksize,chunkmax);
		cchunk = *chunksize-1;
	}
	// This line caused the most pain
	(*c+cchunk)->voxel[chunkPos((vx%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vy%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE, (vz%CHUNKSIZE+CHUNKSIZE)%CHUNKSIZE)] = block;

}

// Loads and unloads world chunks
int loadWorldChunks(int cx, int cy, int cz, Chunk *chunks, int chunkSize, char worldType(int,int,int)) {
	// Unload unneeded chunks
	int shift = 0; // Shift everything this amount
	int chunkRad = CHUNKRAD;
	int chunkMax = (2*chunkRad+1)*(2*chunkRad+1)*(2*chunkRad+1) * 3;
	for(int i=0;i<chunkSize;i++){
		if(chunks[i].x< (cx-chunkRad) || (chunks[i].x> (cx+chunkRad))) { shift++; continue; }
		if(chunks[i].y< (cy-chunkRad) || (chunks[i].y> (cy+chunkRad))) { shift++; continue; }
		if(chunks[i].z< (cz-chunkRad) || (chunks[i].z> (cz+chunkRad))) { shift++; continue; }

		// Shift all of the chunks
		chunks[i-shift] = chunks[i];
	}
	// Load new chunks
	int newSize = chunkSize - shift;
	int newerSize = newSize;
	for(int z=cz-chunkRad;z<=cz+chunkRad;z++) {
		for(int y=cy-chunkRad;y<=cy+chunkRad;y++) {
			for(int x=cx-chunkRad;x<=cx+chunkRad;x++) {
				// Check if chunk already exists
				int chunkExists = 0;
				for(int e=0;e<newSize;e++) {
					if((chunks[e].x==x && chunks[e].y==y) && chunks[e].z==z) {
						chunkExists = 1;
						break;
					}
				}
				if (chunkExists != 0) continue;
				// Load the new chunk

				if(newerSize < chunkMax) {
					chunks[newerSize++] = loadChunkWorld(x,y,z,worldType);
				}
				else {
					exit(-1);
				}
			}
		}
	}
	return newerSize;
}
int main(int argc, char** argv) {
	if(argc != 5 && argc != 6) return -1;

	sscanf(argv[1],"%d", &Width);
	sscanf(argv[2],"%d", &Height);
	sscanf(argv[3],"%d", &MAXSTEP);
	sscanf(argv[4],"%d", &CHUNKRAD);
	if(argc == 5) {
		sscanf(argv[5],"%d", &world_seed);
	}

	// Setup
	if(!(pixels=(unsigned char*)calloc(Width*Height, sizeof(unsigned char)))) {
		printf("Error: Memory Skill issue.\n");
		return -1;
	}

	// Player Stuff
	int use_flashlight = 0;
	float material = 1;
	float speed = 0.1f;

	float px = 2734.0f;
	float py = 8.5f;
	float pz = 850.0f;

	float vx = 0.0f;
	float vy = 0.0f;
	float vz = 0.0f;

	float rx = 0.0f;
	float ry = 0.0f;
	float vrx = 0.0f;
	float vry = 0.0f;

	float height = 1.5f; // Player size
	float width = 0.2f;

	int blockX = 0;
	int blockY = 0;
	int blockZ = 0;
	int bnx = 0;
	int bny = 0;
	int bnz = 0;
	int cx = floorf(px/((float)(CHUNKSIZE)*VOXELSIZE));
	int cy = floorf(py/((float)(CHUNKSIZE)*VOXELSIZE));
	int cz = floorf(pz/((float)(CHUNKSIZE)*VOXELSIZE));
	vec3 sun = normalize((vec3){0.8f, 1.6f, -0.2f});

	char (*daWorld)(int,int,int) = hillyWorld;

	// Chunks
	int cRad = CHUNKRAD;
	int cMax = (2*cRad+1)*(2*cRad+1)*(2*cRad+1) * 3;

	int worldChunkSize = 0;
	Chunk *worldChunks = (Chunk*)malloc(cMax*sizeof(Chunk));
	worldChunkSize = loadWorldChunks(cx,cy,cz,worldChunks,worldChunkSize, daWorld);

	int playerChunkSize = 0;
	int playerChunkMax = cMax;
	Chunk *playerChunks = (Chunk*)malloc(playerChunkMax*sizeof(Chunk));

	struct termios attr;
	tcgetattr(STDIN_FILENO, &attr);
	attr.c_lflag &= ~(ICANON | ECHO);
	tcsetattr(STDIN_FILENO, TCSANOW, &attr);

	unsigned char input[16];
	ssize_t amount;
	struct pollfd pstdin = {.fd=STDIN_FILENO, .events=POLLIN}; // using ".fd" and stuff allows you to define the specific parameter in the struct
	int quit = 0;

	struct timespec start, end;

	// Render
	vec3 looking = (vec3){0.0f,0.0f,0.0f};

	while(quit == 0) {
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start);
		setpos(0,0);
		for(int h=0; h<Height;h++) {
			for(int w=0; w<Width; w++) {
				float u = 2.0f*(float)(w-Width/2)/(float)Height;
				float v = 3.0f*(float)((Height-h)-Height/2)/(float)Height; // flip y

				vec3 ro = {px,py,pz};
				vec3 rd = {u,v,1.0f};
				vec2 rotX = Rot(rd.y,rd.z,rx);
				rd.y = rotX.x;
				rd.z = rotX.y;
				vec2 rotY = Rot(rd.x,rd.z,ry);
				rd.x = rotY.x;
				rd.z = rotY.y;
				vec3 n;
				int vox = 0;

				int bx = 0;
				int by = 0;
				int bz = 0;
				float d = voxelTrace(ro,rd,worldChunks, worldChunkSize, playerChunks, playerChunkSize, &n, &vox,&bx,&by,&bz);

				if(w==Width/2 && h==Height/2) {
					looking = normalize(rd);
					blockX = bx;
					blockY = by;
					blockZ = bz;
					bnx = n.x;
					bny = n.y;
					bnz = n.z;
				}
				vec3 p = vec3_add(ro,vec3_mulf(rd,d));
				// *** Light calculations ***

				// Sun
				float light = getLight(n,sun);
				if(py < -10.0f) {
					light *= 0.6f;
					//light += 3.0f*getLight(n,vec3_mulf(looking,-1.0f))/(d*d);
				}
				if(use_flashlight == 1) {
					light += 3.0f*getLight(n,normalize(rd)) / (d*d);
				}
				light *= 15.0f;
				if(light > 15.0f) light = 15.0f;

				if(d>=MAXDIST) {
					if(dot(normalize(rd),sun) > 0.99f) draw(w,h,11,2); // Sun
					else {
						if(ro.y>-4.0f) draw(w,h,14,2); // Sky
						else draw(w,h,0,0);
					}
				}
				else {
					if((bx==blockX && by == blockY) && bz==blockZ) draw(w,h,7,light);
					else draw(w,h,vox,light);
				}
			}
		}
		draw(Width/2, Height/2,material,15);
		update();
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end);
		double elapsed = (end.tv_sec-start.tv_sec) * 1e6 + (end.tv_nsec-start.tv_nsec) / 1e3;
		setpos(0,Height);
		printf("%f ms   (%d fps)",elapsed*0.001f,(int)(1.0f/(elapsed*0.000001f)));

		// Input
		if(poll(&pstdin, 1, 0) > 0 && (pstdin.revents & POLLIN)) { // If there's data in stdin (0 timeout)
			if((read(STDIN_FILENO, input, 16)) > 0) {
				switch(input[0]) {
					case 9: // TAB
						quit = 1;
						break;
					case 'e':
						material++;
						if(material > 15) material = 1;
						break;
					case 'q':
						material--;
						if(material < 1) material = 15;
						break;
					case 'w':
						vz = speed*cos(ry);
						vx = speed*sin(ry);
						break;
					case 's':
						vz = -speed*cos(ry);
						vx = -speed*sin(ry);
						break;
					case 'd':
						vz = -speed*sin(ry);
						vx = speed*cos(ry);
						break;
					case 'a':
						vz = speed*sin(ry);
						vx = -speed*cos(ry);
						break;
					case 'f':
						use_flashlight = (~use_flashlight)&1;
						break;
					case 32:  // Space
						vy=0.25;
						break;
					/*case 'q':
						vy=-speed;
						break;*/
				case 'c':
						setBlock(material, blockX+bnx, blockY+bny, blockZ+bnz, &playerChunks, &playerChunkSize, &playerChunkMax);
						break;
					case 'x':
						setBlock(AIR, blockX, blockY, blockZ, &playerChunks, &playerChunkSize, &playerChunkMax);
						break;
					case 27:
						if(input[1] == 91) {
							switch(input[2]) {
								case 65:  // Up
									vrx = 0.05f;
									break;
								case 66:  // Down
									vrx = -0.05f;
									break;
								case 67:  // Right
									vry = 0.05f;
									break;
								case 68:  // Left
									vry = -0.05f;
									break;
							}
						}
						break;
				}
			}
		}

		// Update Position/Collision
		float time_cheddar = 0.0001f*elapsed;
		if(time_cheddar > 1.0) time_cheddar = 1.0; // When run onn routers, the elapsed time may be too long and you end up falling through the floor

		if(getBlock(px+vx*time_cheddar+((vx<0.0f)?-width:width),py-height,pz,worldChunks,worldChunkSize,playerChunks,playerChunkSize) == AIR && getBlock(px+vx*time_cheddar+((vx<0.0f)?-width:width),py,pz,worldChunks,worldChunkSize,playerChunks,playerChunkSize) == AIR) px+=vx*time_cheddar;
		if(getBlock(px,py-height,pz+vz*time_cheddar+((vz<0.0f)?-width:width),worldChunks,worldChunkSize,playerChunks,playerChunkSize) == AIR && getBlock(px,py,pz+vz*time_cheddar+((vz<0.0f)?-width:width),worldChunks,worldChunkSize,playerChunks,playerChunkSize) == AIR) pz+=vz*time_cheddar;
		if(getBlock(px,py-height+vy*time_cheddar,pz,worldChunks,worldChunkSize,playerChunks,playerChunkSize) == AIR && getBlock(px,py+vy*time_cheddar,pz,worldChunks,worldChunkSize,playerChunks,playerChunkSize) == AIR) py+=vy*time_cheddar;
		else vy=0.0f;

		if(fabs(rx+vrx*time_cheddar)<(3.14159*0.5f+0.01f)) rx+=vrx*time_cheddar;
		ry += vry*time_cheddar;

		// Add a bit of positional and rotational friction to be more smooth
		vx -= ((vx<0.0f)? (-speed*0.05f*time_cheddar) : (speed*0.05f*time_cheddar));
		vz -= ((vz<0.0f)? (-speed*0.05*time_cheddar) : (speed*0.05f*time_cheddar));
		vy -= 0.02f*time_cheddar;
		vrx -= ((vrx<0.0f)? (-speed*0.03f*time_cheddar) : (speed*0.03f*time_cheddar));
		vry -= ((vry<0.0f)? (-speed*0.03f*time_cheddar) : (speed*0.03f*time_cheddar));

		if(fabs(vx) < speed*0.2) vx = 0.0f;
		if(fabs(vz) < speed*0.2) vz = 0.0f;
		if(fabs(vrx) < speed*0.1) vrx = 0.0f;
		if(fabs(vry) < speed*0.1) vry = 0.0f;
	}
	return 0;
}
