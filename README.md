# Generator for 3D Reuleaux/Constant-Width Shapes
Generate a point cloud of points on the surface a random reuleaux/constant width shape
___
## Warning
No compiled version yet. To run you first need to download the [eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page) to the 3rdparty folder. The final path should be: "/3rdparty/eigen-eigen-xxxxx/Eigen".
Then open ConstantWidthShapes.pro using [Qt](https://www.qt.io/).
Also, this project startet just as a first execise in Qt, but it turned out to actually produce working models. So dont expect great code to begin with.

## From Output to 3D mesh
The output file will be a .obj file containing only vertices. To turn this in to a 3D mesh one can use [MeshLab](http://www.meshlab.net/).
First import the OBJ file then use: Filters>Remeshing, Simplification and Reconstruction>Convex Hull. Then you can export the mesh in your desiered fomat.

9 example meshes generated using this method can be found here: https://drive.google.com/file/d/1iWjz1TvAEv5ugF3FaL8wvvbQ-em-y5bG/view?usp=sharing

![alt text][ImageOfExample]

[ImageOfExample]: https://github.com/oyboy/Generator_for_3D_Reuleaux_or_Constant_Width_Shapes/blob/master/imgs/exampleShapes.png "The 9 example meshes"
A simulation with these examples: https://youtu.be/qGcFCXV-Vsg
## How?
We start with only working with paris of points where the distance between the points in each pair is always the same.


- The first step is generating *n* such pairs with a random position and rotation
- We then iterate the position of each pair untill it satisfies the following condition then add it to a set of valid pairs.

For a pair to be valid: From each point in the pair, no point in the set of already valid pairs can be further away than the other point in the pair. In other words when we add a pair to the set, form each point in the pair the point in the set furthest away should be the other point in the pair.

This might seem like a simple rule but the fact that it seems to always generate constant width shapes raises some questions that I dont find imediatly obvious:
- Why would this limitation generate convex hulls? It makes sence in a local area but also the object as a whole is convex.
- These shapes are only of constant width if the two points in contact, when placed betwen two parallel planes, form a line perpendicular to the planes. Will this always be the case for all constant width shapes?
- If yes, would this be able to generate all constant with shapes(in 3 dimensions)?

## Improvements
The iteration could be replaced by a more dererministic approach if one can show that a pair will be placed near pairs with a similar angle(as far as I have seen this has always been the case for the shapes generated). One might then be able to find all possible positions of the pair by just looking at the three clossest pairs(in angle) that enclose the pair we want to add. Currently the iteration process just uses the first valid position it reaches.

Imrove the tranformation between constant width shapes so that one can specify what shape to converge to. The current result is random as seen below.

![alt text][AnimExample]

[AnimExample]: https://github.com/oyboy/Generator_for_3D_Reuleaux_or_Constant_Width_Shapes/blob/master/imgs/Anim1.gif "Animation Example"
## Related Papers
Here is a list of papers related to the topic.
None of the methods showcased here are similar to this one as far as I have understood, though the first paper produces results close(at least on a subjective level, just looking at the shapes) to mine.

[MEISSNER POLYHEDRA](https://arxiv.org/pdf/1608.06354.pdf)  
LUIS MONTEJANO, EDGARDO ROLDAN-PENSADO

[Meissner’s Mysterious Bodies](https://www.researchgate.net/publication/225748121_Meissner's_Mysterious_Bodies)  
Bernd Kawohl, Christof Weber

[Bodies of constant width in arbitrary dimension](https://hal.archives-ouvertes.fr/hal-00385113/document)  
Thomas Lachand-Robert, Edouard Oudet


[Parametric Shape Optimization using the Support Function](https://arxiv.org/pdf/1809.00254.pdf)  
Pedro R.S. Antunes, Beniamin Bogosel


[Shape Optimization Under Width Constraint](https://link.springer.com/article/10.1007/s00454-012-9471-z)  
Édouard Oudet


[MAXIMAL R-DIAMETER SETS AND SOLIDS OF CONSTANT WIDTH](https://arxiv.org/pdf/1003.5824.pdf)  
ETHAN AKIN


[On Curves and Surfaces of Constant Width](https://arxiv.org/pdf/1504.06733.pdf)  
H. L. Resnikof

Article proving the existance of a contious transformation between constant width shapes exist so that every step in the transformation is also of constant width.
https://www.researchgate.net/publication/225624683_Topology_of_the_hyperspace_of_convex_bodies_of_constant_width  
I would not claim to understand the proof.
