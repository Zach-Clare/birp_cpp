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
- `refine/` branches will reduce or eliminate rendering artefacts, thereby "refining" the image.
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

### Generating rays in camera-space
Okay, how does this work? The main idea revolves around an image plane placed 1 unit away from the camera. The program is given an FOV from which it can calculate how tall (y) and wide (x) this image plane is. It is also provided with an image dimension in pixels. The centre pixel (if it exists, which it won't if there are an even number of pixels) will be pixel 0, 0. The x value increases to the right, and the y value increases upwards. So, the top left pixel has the coordinate -1, 1, as it is left (-1) and up (1). The top right is 1, 1 because it is positive in both directions. The bottom left is -1, -1 and finally the bottom right is 1, -1. Remember that this plane is physically placed 1 unit away from the camera in 3d space. By convention, the plane is located 1 unit along the -z (negative z) axis. It is possible to place it 1 unit away in the positive z, or the negative x, or anywhere for that matter. You just have to account for it elsewhere. In this program, like many many others, we will be placing it -1z away.

Each pixel needs a vector (called a ray in computer graphics, which is what we will be calling it from now on) that originates from the camera and goes through the centre of that specific pixel on the plane. This means that we will be generating H x W amount of rays. To generate this vector, we need to know the origin of the camera and the location of that pixel in space. We are defining a local space around the camera, so the camera is 0, 0, 0. We can get the x, y coordinates of the pixel by first normalising the pixel number into [0, 1] space. To do this, we add 0.5 to the pixel index (ensuring we're calculating the position of the centre of the pixel, not where the pixel begins), and divide it by our total number of pixels along that dimension. So for the fourth pixel along our top row in an image that's 144 pixels wide (remembering that arrays start at zero), we do (3 + 0.5) / 144 = 0.024305 which is our x value. Since the pixel is in the top row (which is our final row if this coordinate system begins in the bottom left), we can do (143 + 0.5) / 144 = 0.996527. To convert this into our [-1, -1] space, we can simply double it and minus 1. That puts those coordinates at -0.951388 and 0.993055 respectively.

We could use this to get our pixel ray, but we need to account for the FOV first. Notice that the FOV is the angle between the each side of the image (usually referring to the right and left, but since our image here is square, it refers to the angle between each side of the image and it's opposite side). Notice, too, that half of that FOV angle will take us from the middle to one of the sides. We've touched already about how the centre pixel's coordinates will be 0, 0, -1. I hope you can see how the distance covered by that vector is simply 1 unit. The plane is at a right angle to this centre vector and goes upwards (and downwards) 1 unit. We can construct an ad-hoc right-angle triangle with the line between the camera and the plane acting as our adjacent side and the top holf of the image plane acting as our opposite side. Here, the pixel ray will act as teh hypotenuse. Thankfully, because our opposite and adjacent sides are 1 unit, it's as simple as multiplying the coordinate value with _tan(FOV / 2)_.

If any of this didn't make sense, or you'd like a more in-depth explaintion, I found [ScratchPixel's explaination for camera ray generation](https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays.html) really useful.

So that's great, that's fantastic, we have a ray that describes where to go from the camera! But there's a bit more to go before we can start integrating along this ray.

### Converting those rays to world-space
These rays are accurate in relation to the camera, but unfortunately we use camers to look at things. This means we need to put it in the right place and point it where we want to see. The placement of the camera is a fairly simple translation, but this pointing business is a bit more tricky. I originally used a spherical rotation system and converted between spherical and cartesian, but the rays were squashed together toward the bottom and the top of the view (as spherical coordinates bulge at the equator). I then came up with another system of approaching the rotation of the camera as two seperate 2D rotations. I really wanted that to work and I drew a _lot_ of diagrams to wrap my head around it, but in the end I figured it was best to use something that was industry-standard. There would be support for it and, with any luck, programmers and physicists (and anyone else for that matter) that read the code might be able to recognise it from elsewhere and that may aid their undertanding.

The method used in ths program is to calculate the rotation angles, plug them into a set of well-defined rotation matrices, unify those matrices into one rotation matrix, and perfrom matrix multiplication on our camera ray with that rotation matrix to find the new ray direction. The best way I've come to understand what our rotation angles are specifically is - well let's say our camera is pointed exactly where we want it at the very beginning, how do we rotate the camera so that it's orientation agrees with GSEs definition of where x, y, and z are? I realise there's a bit going on here so let's go one step at a time.

Let's go through an example here - we'll use a specific moment in orbit where the position of the spacecraft in xyz GSE coordinates are `[5.9, 6.7, 17.7]`. For this orbit position, we are given the aimpoint of `[7.8, 0, 0]`. So for perspective, the spacecraft is quite high up above the earth, and is in front and to the side. We are finding the angles to rotate along each axis if we reverse the camera's orientation when the negative z axis points directly at the aimpoint which lies at some point between the sun and the earth. 

Let's define our first angle of rotation along the spacecraft's x axis. Remember that the spacecraft's x axis goes from left to right across the image plane. Rotating along the x axis looks (from the camera's perspective) as if we are panning the shot up and down. Here is a a not-to-scale diagram showing the value we're trying to get. By using the coordinates available to us already, we can create another oh-so-handy ad hoc right-angle triangle (shown in green). By constructing this right-angle triangle and using SOHCAHTOA rules, we can find the angle we're interested in. We know the adjacent side, as that is simply the spacecraft's z coordinate. In order to find the opposite length (along the bottom), we'll need to create another triangle. This new triangle is a right-angled triangle along the x-y plane flat as if on the floor (shown in blue). Notice that the opposite length we're trying to find here is actually the hypotenuse of our newer on-the-floor-triangle. By squaring, summing, and then square root-ing these lengths, we can find the legnth of that hypotenuse and, therefore, the opposite length of our stood-up triangle with arctan.

![3noscale_enhanced_colour](https://github.com/Zach-Clare/birp_cpp/assets/41343750/ad73bb19-0aeb-4fac-96f0-87956ce0cde9)

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
Now we have the opposite length, we use $tan$ from SOHCAHTOA on the $opposite$ over the $adjacent$:
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

There we have our $R_{x}(\theta)$. Now that we've defined a rotation matrix that rotates our vectors and rays around their x axis, I want to explain where I went wrong on the next step where we calculate rotation along the y axis.

Now imagine that we are floating directly above the earth looking down through GSEs negative z axis. We need to align the x and y axis of the spacecraft to the GSE x and y axis, but the y axis is currently pointing at the aimpoint still and the x axis is 90 degrees right of that. Remember that we are trying to reverse this transformation so it goes back to GSE. To get y back to GSEs, see that we can define that with the same figures we used to calculate the hypotenuse of the flat blue traingle. Namely:
```math
\begin{aligned}
&&R_{z}(\theta) = &\arctan{\frac{1.9}{6.7}}\\ 
&&              = &0.27632 radians\\ 
&&              = &15.8323 degrees\\ 
\end{aligned}
```

Looking at this, it's easy to see why that seems reasonable, but if we use these values in the rotation matrix, we'll see we slightly overcook it and pan too far up the x axis. If we take into account the height of the spacecraft during this calculation, you can see how the angle might become significantly more acute. Here is another not-to-scale diagram:

![smile orbit rotations](https://github.com/Zach-Clare/birp_cpp/assets/41343750/029aa458-02bc-47e8-a3dd-8f4a18f2bd95)

Our previous on-the-floor triangle (in a light blue) has a wider angle under the spacecraft than our new orange triangle. How do we calculate the angle at the orange triangle? It might be a bit hard to notice from this angle (sorry) but that orange triangle is a right-angle triangle (again!). By finding the length along the side and the distance along the x axis (the latter of which we already have) we can find the smaller angle that we need. To find the length of the adjacent side, we can use pythagoras theorum and we just square, sum, and square root the lengths. That looks like this:
```math
\begin{aligned}
&&hypotenuse = &\sqrt{spacecraft_ {y} ^2 + spacecraft_ {z}^2}\\
\\
&&            = &\sqrt{6.7^2 + 17.7^2}\\
\\
&&            = &\sqrt{44.89 + 313.29}\\
\\
&&            = &\sqrt{358.18}\\
\\
&&            = &18.92
\end{aligned}
```

And then we use that in yet another $\tan$ calculation:
```math
\begin{aligned}
&&\theta = &\arctan{\frac{opposite}{hypotenuse}}\\
\\
&&       = &\arctan{\frac{1.9}{18.92}}\\
\\
&&       = \arctan{0.100422}\\
\\
&&       = 0.10008728 radians\\
\\
&&       = 5.73457881 degrees\\
\end{aligned}
```
And if we use that angle instead of our earlier mis-calculated one, we'll end up with a perfect setup for another rotation matrix. On that note, he's the rotation matrix for a rotation about the z axis:


---------------

By defining the rotations with a series of rotational matrices, we can calculate a single matrix and perform the rotation for our ray in a single step. The [wikipedia page for 3D rotations](https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions) contains a really useful record of the rotational matrices for each axis of rotation, and by calculating the angles to rotate by in each axis (or maybe only a couple), we can calculate the rotation matrix specific for this position in orbit and aimpoint. Then, through the magic of matrix multiplication (with caveats), we can calculate the new vector for each ray.

We will label the line connecting the spacecraft's x-y projection (the point below the spacecraft) to the main x axis (whose length is equivalent to the spacecraft's y coordinate) the adjacent side. That would make the side where y = 0 the opposite side. 
