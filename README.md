# Terminal Based Voxels
A Terminal-based Voxel Raycaster that will run on literally anything (including a router)

## The Goal
The goal of this project was to make an actual text-based Minecraft clone using only the Linux terminal
and the C standard library. This way it would be hopefully be able to be compiled and ran on basically any (Linux) system.

## Technical Parts
Most of the technical parts of this project took inspiration from the GNU nano text editor and "top" command.
The use of nano keybinds showed that it was possible to get keyboard input without hitting the enter key and the
"top" command showed that it was possible to poll keyboard input without having the whole program wait for the 
user to press something.

In the end, the answer was simply a special stdin mode that stops it from stalling the program and from needing to press
the enter key.

## Compiling
For the minimal Voxel Game:
gcc voxel_min.c -lm -o voxel_min

## Usage
For the minimal Voxel Game:
run with ./voxel_min [screen width (in text characters)] [height (in text characters)] [max render distance (100 is usually fine)] [chunk radius (chunks to load in all directions)] [optional: world seed]

WASD to move around
C to place the current block in front of you (indicated by cursor in middle of screen)
X to break the block in front of you
E and Q to switch current block
SPACE to jump

Note that you can only press one key at a time as it still uses stdin but the movements have some deacceleration time to account for that

## Problems running on a router
There were many issues with getting this to run on a router but it was mainly due to three problems:
1. The router was an old MIPS CPU with no FPU (no floats or sqrt function)
2. There was no available cross-compiler for the CPU's architecture
3. The router took way too long to do anything and needed to be power cycled every time because of it 
