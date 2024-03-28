# SAC-for-H-Bond-Learning
Code currently adopted from https://github.com/p-christ/Deep-Reinforcement-Learning-Algorithms-with-PyTorch

Paper: https://arxiv.org/pdf/2010.02193.pdf

Current SAC Discrete Test on Cart Pole Environment

![alt text](https://user-images.githubusercontent.com/47870060/116022149-78174800-a5fe-11eb-995d-fece01c3e172.png)

# Running Octree
To run the octree, go to `/Octree/` and run `make octree`. Then run `python Octree.py` to run the current version of the Octree

# TODOs
Current the entire Octree is wrapped but the individual functions are not. Specifically, the process of removing atoms and inserting atoms is not defined in the `cython` code.