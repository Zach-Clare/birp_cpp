-----------------------------------------------------------
    This file contains the Unpublished Intellectual Property of
    University College London and All Rights Are Reserved.
    Copyright (c) University College London, 2024
-----------------------------------------------------------


# The Boundary Image Rendering Project
Welcome to the Boundary Image Rendering Project (BIRP). This project simulates clean and unprocessed renderings of a simulated magnetopause as part of the data pipeline for the Solar wind Magnetosphere Ionosphere Line Explorer (SMILE) mission. These simulated images are what the Soft X-Ray Imager (SXI) intrument would see before any image degredation. The images generated by this project will be sent through various intrument simulators and then to a machine learning model so that, when the mission launches in early-mid 2025, accurate estimations can be made about our magnetosphere's location in real time.

This project was written for University College London (UCL) by Zach Clare. It is based on an IDL implementation authoured by Dr. Andrew Read and Dr. Steven Sembay at the University of Leicester (UoL).

The SMILE project is a joint [European Space Agency (ESA)](https://www.esa.int/) and [Chinese Academy of Sciences (CAS)](https://english.cas.cn/) mission.

**In laymen's terms**, this software simulates an X-Ray camera orbiting the Earth pointing at our magnetic field. This helps plasma scientists study how our magnetic field is influenced by the sun's solar wind and forms part of a space research mission launching in December 2025.

## Code Structure

See `main.cpp` for a good example on using this software. If you don't fancy doing that, all you need to know is that there are two main components - `DataCube` and `Camera`. `DataCube` describes the world in which the simulation is ran in and holds the data from the Magneto-Hydro Dynamic (MHD) simulation. `Camera` holds all the logic for positioning and rendering. `Helper` (named `Handbook` in the Python version of this project) is a static class that holds a few helper functions that don't quite belong in the main two classes.

## Repository Organisation

Just some quick organisational points that keep us all on the same page.

### Branches

The `main` branch should always hold the most stable and reliable version of BIRP. Branch names may be prefixed with a keyword:
- `speed/` branches hold new optimisation methods that will "speed" up the rendering process, for example `speed/matrix_multiplication` might describe a branch that enhances the speed certain calculations by utilising matrix operations.
- `refine/` branches will reduce or eliminate rendering artefacts or improve scientific accuracy, thereby "refining" the image.
- `feature/` branches introduce a new non-science "feature" (or features). For example, `feature/position_as_day` could describe a branch that adds support for specifying the position of the camera as a day of orbit rather than an [x,y,z] coordinate, and `feature/HDF5` might add support for exporting as an HDF5 image.
- `fix/` branches "fix" bugs that aren't necessarily science-related, though they may be related to science components of the code.

### Development Flow

When creating a branch, make sure the `main` branch is checked out first. This will use the `main` branch as a basis for the new branch. If you are wanting to explore changes to an existing child branch, check out that child branch first to create the new branch from there.

To implement changes from a development branch into `main`, keep merging the development branch back into it's parent branch until you reach main.

For example, the typical development flow may look like this:
- We want to implement tricubic interpolation for more acurate renders.
- So we begin from the `main` branch with `git checkout main`.
- A development branch can be created with `git checkout -b refine/tricubic`.
- This is where we would write our code, perform tests, make sure it's all working.
- If we want this to become part of the core simulator, we need to merge it back up the tree.
- Make sure all your changes are pushed to the development branch with `git status`. If there are outstanding changes, commit and push with an appropriate commit message.
- Make sure we have the latest version of main with `git checkout main` and then `git pull`.
- If there are new changes, update the development branch to ensure your changes are still working properly with `git checkout refine/tricubic` and `git merge main`.
- Perform your final tests.
- If all is well, `git checkout main` and `git merge refine/tricubic` will introduce your changes. `git push` will push your changes to the repository for all to see.
- Sit back an celebrate your genius.

## The method
This area describes (without code) the maths behind how this program works. There's a bit of trigonometry, but nothing too intense. With that all said, I'm a programmer - not a mathemetician. If you see a failing in this method please get in touch and discuss it with me. However, this method does produce demonstratably good results and I'm confident this is correct.

> **Images and mathematics may not display correctly in the Doxygen documentation. For best results, read these docs through the [GitHub repo](www.github.com/Zach-Clare/birp_cpp).**

### Generating rays in camera-space
Okay, how does this part work in a nutshell? The main idea revolves around an image plane placed 1 unit away from the camera. The program is given an FOV from which it can calculate how tall (y) and wide (x) this image plane is. It is also provided with an image dimension in pixels. The centre pixel (if it exists, which it won't if there are an even number of pixels) will be pixel 0, 0. The x value increases to the right, and the y value increases upwards. So, the top left pixel has the coordinate -1, 1, as it is left (-1) and up (1). The top right is 1, 1 because it is positive in both directions. The bottom left is -1, -1 and finally the bottom right is 1, -1. Remember that this plane is physically placed 1 unit away from the camera in 3d space. By convention, the plane is located 1 unit along the -z (negative z) axis. It is possible to place it 1 unit away in the positive z, or the negative x, or anywhere for that matter. You just have to account for it elsewhere. In this program, like many many others, we will be placing it -1z away.

![grids](https://github.com/user-attachments/assets/91d3db49-aefc-464e-92c9-b8842fcc14f9)

To explain more clearly, each pixel needs a vector (called a ray in computer graphics) that originates from the camera and goes through the centre of that specific pixel on the image plane. This means that we will be generating H x W amount of rays. To find out where on the -1 to +1 scale our pixel is, we need to do a few things. We first need to know the origin of the camera - we are defining a local space where the camera is the origin, so the camera is 0, 0, 0. We can get the x, y coordinates of the pixel by first normalising the pixel number (as in, pixel number 1, number 2, etc. until we finish that row) into [0, 1] space. To do this, we add 0.5 to the pixel index (ensuring we're calculating the position of the centre of the pixel, not where the pixel begins), and divide it by our total number of pixels along that dimension. So for the fourth pixel along our top row in an image that's 144 pixels wide (remembering that arrays start at zero), we do (3 + 0.5) / 144 = 0.024305 which is our x value. Since the pixel is in the top row (which is our final row if this coordinate system begins in the bottom left), we can do (143 + 0.5) / 144 = 0.996527. To convert this into our [-1, -1] space, we can simply double it and minus 1. That puts those coordinates at x-0.951388 and y0.993055.

We could use this to get our pixel ray, but we need to account for the FOV first. Notice that the FOV is the angle between the each side of the image (usually referring to the right and left, but since our image here is square, it refers to the angle between each side of the image and it's opposite side). Notice, too, that half of that FOV angle will take us from the middle to one of the sides. We've touched already about how the centre pixel's coordinates will be 0, 0, -1. One can see how the distance covered by the vector describing the camera origin to the centre pixel is simply 1 unit. The plane is at a right angle to this centre vector and goes upwards (and downwards) 1 unit. We can construct an ad-hoc right-angle triangle with the line between the camera and the plane acting as our adjacent side and the top half of the image plane acting as our opposite side. There's a great diagram linked just below. Here, the pixel ray will act as the hypotenuse. Thankfully, because our opposite and adjacent sides are 1 unit, it's as simple as multiplying the coordinate value with _tan(FOV / 2)_.

Once each coordinate (apart from the z coordinate) is multiplied by this factor, we're left with a 3D vector. That vector, originating from the origin point of our camera (which is not the same origin as our pixel grid, which would be the earth), describes the line of sight for each pixel of our grid. That's great, but leads us onto our next problem...

Just as a quick note before we move on, if any of this didn't make sense, or you'd like a more in-depth explaintion, I found [ScratchPixel's explaination for camera ray generation](https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays.html) really useful. The diagrams especially are very useful for understanding this concept.

### Converting those rays to world-space
These rays are accurate in relation to the camera, but unfortunately we need to put the camera in the right place and point it where we want to see. The placement of the camera is a fairly simple translation, but this pointing business is a bit more tricky. I originally used a spherical rotation system and converted between spherical and cartesian, but the rays were squashed together toward the bottom and the top of the view (as spherical coordinates bulge at the equator). I then came up with another system of approaching the rotation of the camera as two seperate 2D rotations. I really wanted that to work and I drew a _lot_ of diagrams to wrap my head around it (and I've edited this readme.md a lot of times because of that), but in the end it's working with rotation matricies and two 2D rotations.

In a nutshell, the method used in ths program is to calculate the rotation angles, plug them into a set of well-defined rotation matrices, unify those matrices into one rotation matrix, find the inverse of that matrix, and perfrom matrix multiplication on our camera ray with that rotation matrix to find the new ray direction. The best way I've come to understand what our rotation angles are specifically is - well let's say our camera is pointed exactly where we want it to be pointed at the very beginning, how do we rotate the camera so that it's orientation agrees with GSEs definition of where x, y, and z are? I realise there's a bit going on here so let's go one step at a time.

We can use an example - we'll use a specific moment in orbit where the position of the spacecraft in xyz GSE coordinates are `[5.9, 6.7, 17.7]`. For this orbit position, we are given the aimpoint of `[7.8, 0, 0]`. So for perspective, the spacecraft is quite high up above the earth, and is in front and to the side. We want to find how much we should rotate along each of the spacecraft's axis so that both sets of axes agree. 

![correct_rotations_colour](https://github.com/user-attachments/assets/b3b15c88-12e2-43e0-896a-2c63f1c8b8e5)

Let's focus on the diagram above and define our first angle of rotation along the spacecraft's x axis. Remember that the spacecraft's x axis goes from side to side across the image plane. Rotating along the x axis looks (from the camera's perspective) as if we are panning the shot up and down. By using the coordinates available to us already, we can create a right-angle triangle (shown in blue) that describes the angle of rotation between the spacecrafts x axis and GSE's x axis. By constructing this right-angle triangle and using SOHCAHTOA rules, we can find the actual angle. We know the adjacent side, as that is simply the spacecraft's z coordinate - it's height. In order to find the opposite length (along the bottom, in the XY plane), we'll need to create another triangle. This new triangle is a right-angled triangle along the XY plane flat as if on the floor (shown in pink). Notice that the blue triangle's opposite length we're trying to find here is actually the hypotenuse of our newer pink triangle. By squaring, summing, and then square root-ing these lengths, we can find the legnth of that hypotenuse and, therefore, the opposite length of our stood-up triangle with arctan. 

$$opposite = \sqrt{(x_ {aim} - x_ {spacecraft})^2 + (y_ {aim} - y_ {spacecraft}) ^2}$$

$$\angle z_{s}z = \arctan({\frac{opposite}{adjacent}})$$

That's great, we now have our angle to rotate around x. If we refer to [wikipedia page for 3D rotations](https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions), we can see that the rotation matrix for any 3D rotation around the x axis (ever!) looks like this:

```math 
R_{x}(\theta) = \begin{bmatrix}
       1 &           0 &             0  \\[0.3em]
       0 & \cos{\theta} & -\sin{\theta} \\[0.3em]
       0 & \sin{\theta} & \cos{\theta}
     \end{bmatrix}
```
If we were to calculate this for our example, where $\theta$ is $\angle z_{s}z$:
```math
\begin{aligned}
&&opposite = &\sqrt{(7.8 - 5.9)^2 + (0 - 6.7) ^2}\\ 
&& = & \sqrt{(1.9^2 + 6.7^2}\\ 
&& = & \sqrt{(3.61 + 44.89}\\ 
&& = & \sqrt{(48.5}\\ 
&& = & 6.96419...\\
\end{aligned}
```
Now we have the opposite length, we use $arctan$ from SOHCAHTOA on the $opposite$ over the $adjacent$:
```math
\begin{aligned}
&&\theta = &\arctan({\frac{6.96419}{17.7}})\\
\\
&& = &\arctan{0.393457}\\
\\
&& = &0.37485342\\
\end{aligned}
```
There we have our angle in radians. To calculate the rotation matrix, simply plug it in and calculate:
```math
\begin{aligned}
&&R_{x}(\theta) = &\begin{bmatrix}
       1 &            0 &             0  \\[0.3em]
       0 & \cos{\theta} & -\sin{\theta} \\[0.3em]
       0 & \sin{\theta} &  \cos{\theta}
     \end{bmatrix}\\
\\
&& = &\begin{bmatrix}
       1 &                0 &                 0  \\[0.3em]
       0 & \cos{0.37485342} & -\sin{0.37485342} \\[0.3em]
       0 & \sin{0.37485342} &  \cos{0.37485342}
     \end{bmatrix}\\
\\
&& = &\begin{bmatrix}
       1 &          0 &           0  \\[0.3em]
       0 & 0.92358911 & -0.38338382 \\[0.3em]
       0 & 0.38338382 &  0.92358911
     \end{bmatrix}
\end{aligned}
```

There we have our $R_{x}(\theta)$. Now that we've defined a rotation matrix that rotates our vectors and rays around their x axis, Let's have a look at the next angle.

Now imagine that we are floating directly above the earth looking down through GSEs negative z axis, after having already reversed our camera's rotation about its x axis. That means the camera'z z now agrees with GSEs z and both of these axes are pointing up out the page towards us. We next need to align the x and y axis of the spacecraft to the GSE x and y axis, but the y axis is currently pointing in the direction of the aimpoint and the x axis is 90 degrees right of that. We'll need to calcuate the rotation through the movement similar to a spinning top. To get the spacecraft's y back to GSE's y, see that we can define that with the same figures we used to calculate the hypotenuse of the flat pink traingle. Namely:
```math
\begin{aligned}
&&R_{z}(\theta) = &\arctan{\frac{1.9}{6.7}}\\ 
&&              = &0.27632 radians\\ 
&&              = &15.8323 degrees\\ 
\end{aligned}
```
And if we use that angle in the same way using the appropriate 3D rotation matrix, we get the below:

```math
\begin{aligned}
&&R_{z}(\theta) = &\begin{bmatrix}
       \cos{\theta} & -\sin{\theta} &  0  \\[0.3em]
       \sin{\theta} &  \cos{\theta} &  0  \\[0.3em]
                  0 &             0 &  1
     \end{bmatrix}\\
\\
&& = &\begin{bmatrix}
        \cos{0.27632} & -\sin{0.27632} &  0  \\[0.3em]
        \sin{0.27632} &  \cos{0.27632} &  0  \\[0.3em]
                    0 &           0 &  1
     \end{bmatrix}\\
\\
&& = &\begin{bmatrix}
       0.96206592 & -0.2728171 &  0  \\[0.3em]
        0.2728171 & 0.96206592 &  0  \\[0.3em]
                0 &          0 &  1
     \end{bmatrix}
\end{aligned}
```

```math

\begin{bmatrix}
       0.962 & -0.273 &  0.00  \\[0.3em]
       0.256 &  0.889 & -0.383  \\[0.3em]
       0.105 & -0.369 &  0.924
\end{bmatrix}

```

Now that we have both sets of rotation matraces, we need to use matrix multiplication to create a single matrix. Now, as I stated earlier, I'm not a mathematician, but I am confident this is correct. As you might know, matraces aren't commutative, meaning it matters which way round you multiply them so we need to pay attention to which one we calculated first in our example. Once you've done the multiplication, you will need to use the inverse of that vector to multiply it by the pixel vector when you want to find that vector in world space, in GSE. The program does $R_{x} \times R_{z}$, finds the inverse and multiplies it with the vector for each pixel, and that seems to work perfectly.

Once we have that resultant vector, we make it into a unit vector and then multiply it by a uniformly increasing set of factors in order to get multiple vectors with specific distances. These are vectors from our spcaecraft to our sample points. We minus the spacecrafts position from each vector which will give us their 3D coordinate. We take the sample at each of these coordinates and add them all to a running total flux level for that pixel. Once we've taken samples from all the vectors, we've reached the end of the ray and we move onto the next pixel. Once we've completed the last pixel, we're done! The SXI unit rotates itself so that the sun-earth line (our x axis) is perpendicular to the sides of the image. There's a great visualisation of this that makes that really easy to understand, but I'm not sure where it is. Once I do, I'll post a link here. To do this, we project the position of the sun onto a 2D plane of the image, and then calculate how much to rotate the image by based on that. That's all in Orient and Point2Plane and needs refactoring for clarity.

## Datacube and CMEM
When passed an `-i` flag, the renderer will step through and sample the datacube passed. When given the `-c` flag, the renderer will use CMEM to calculate emission levels. Below is a full list of CMEM input parameters:
- x, y, z, a – GSE position and aimpoint where aimpoint = (a, 0, 0) 
- s – pixel size in degrees 
- v – solar wind velocity 
- b – magnetic parameter (I may be mistaken here) 
- d – dipole tilt 
- p0, p1, p2, p3 – CMEM parameters 
- q – CMEM B 
- r – CMEM alpha 
- f – CMEM beta 
- e – CMEM bs 
- g – CMEM A1 
- u – CMEM A2 
- j – CMEM bs_ay bowshock flaring 
- k – CMEM bs_az bowshock flaring 

## Reading the data
To see the data in the way we usually want to, we need to convert it to a FITS image. ~~I use a Python program I created, reading from an intermediary data file to get from C++ to Python. I did try creating it as one unit in C++ but it seemed very finnicky.~~ I have since fixed this and included the EleFits library so that FITS exports can all be done in C++. 

# Build Instructions
## Pre-requisites
To run this, you need to have installed <a href="https://github.com/CNES/EleFits">EleFits</a> and <a href="https://heasarc.gsfc.nasa.gov/fitsio/">cfitsio</a>. The CMake file in this project (the file which tells the computer how to compile the code) will look for your installed version of EleFits, and EleFits (I think) will need to look fro your install of cfitsio. It's a bit finnicky so it's probably best you have a bit of knowledge on building C++ projects. Both EleFits and cfitsio have their own instructions for installing their software in the links above.

## Download and build
1. Once EleFits is installed, download this project code code with `git clone` into wherever you want to have the files on your computer.
2. Once downloaded, navigate to wherever you donwloaded the project in your terminal. In the root of the project, create a build directory with `mkdir build`. The directory can be called whatever you want, it's just convention to call it "build".
3. Next we'll prepare the files for compilation. Do this with `cmake -DCMAKE_MODULE_PATH=../EleFits/cmake/modules ..`. We're calling cmake and giving it a flag, telling it where it can find the modules needed. This path, in this project specifically, points to a cmake file that tells cmake where to look for EleFits. It's not an ideal way of doing this but this is down to the EleFits developers as far as I'm aware. If you know a bit about C++ compilation with C++ and know a better way of doing this, please reach out.
4. Now we will actually build the files with `cmake --build .`. Passing the build flag to cmake means we also need to give it a location to put the files which, since we're currently in the build directory, is `.`.
5. That will create a `birp` executable. You can place the batch files in the `data/batch` folder and run `./birp` to run it. Your output FITS files will be in `data/output`.



---------------

By defining the rotations with a series of rotational matrices, we can calculate a single matrix and perform the rotation for our ray in a single step. The [wikipedia page for 3D rotations](https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions) contains a really useful record of the rotational matrices for each axis of rotation, and by calculating the angles to rotate by in each axis (or maybe only a couple), we can calculate the rotation matrix specific for this position in orbit and aimpoint. Then, through the magic of matrix multiplication (with caveats), we can calculate the new vector for each ray.

We will label the line connecting the spacecraft's x-y projection (the point below the spacecraft) to the main x axis (whose length is equivalent to the spacecraft's y coordinate) the adjacent side. That would make the side where y = 0 the opposite side. 
