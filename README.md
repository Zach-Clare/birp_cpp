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
