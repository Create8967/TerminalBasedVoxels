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

## Problems running on a router
There were many issues with getting this to run on a router but it was mainly due to three problems:
1. The router was an old MIPS CPU with no FPU (no floats or sqrt function)
2. There was no available cross-compiler for the CPU's architecture
3. The router took way too long to do anything and needed to be power cycled every time because of it 
