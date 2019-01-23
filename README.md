# Generator for 3D Reuleaux/Constant-Width Shapes
Generate a point cloud of points on the surface a random reuleaux/constant width shape
___
## Warning
No compiled version yet. To run you first need to download the [eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page) to the 3rdparty folder. The final path should be: "/3rdparty/eigen-eigen-xxxxx/Eigen".
Then open ConstantWidthShapes.pro using [Qt](https://www.qt.io/).
Also, this project startet just as a first execise in Qt, but it turned out to accually produce working models. So dont expect great code to begin with.

## From Output to 3D mesh
The output file will be a .obj file containing only vertices. To turn this in to a 3D mesh one can use [MeshLab](http://www.meshlab.net/).
First import the OBJ file then use: Filters>Remeshing, Simplification and Reconstruction>Convex Hull. Then you can export the mesh in your desiered fomat.

9 example meshes generated using this method can be found here: https://drive.google.com/file/d/1iWjz1TvAEv5ugF3FaL8wvvbQ-em-y5bG/view?usp=sharing

![alt text][ImageOfExample]

[ImageOfExample]: https://github.com/oyboy/Generator_for_3D_Reuleaux_or_Constant_Width_Shapes/blob/master/imgs/exampleShapes.png "The 9 example meshes"
A simulation of these examples: https://youtu.be/qGcFCXV-Vsg
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
The iteration could be replaced by a more dererministic approach if one can show that a pair will be placed near pairs with a similar angle. One might then be able to fint all possible positions of the pair by just looking atthe three clossest pairs(in angle) that enclose the pair we want to add.

Add posibility to re-iterate after moving one+ pairs and potentialy making an animation of countiously changing constant width shapes.
