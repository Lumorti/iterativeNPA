
# NPA test suite

This repo contains code I use for testing NPA related things.

## Usage

If you want to test this code (because maybe you're reviewing a paper in which I reference this), you need to install MOSEK, Optim, Eigen and SCS. See the various websites for how to do that, you also need a license (academic or otherwise) for MOSEK.

You then need to set the environment variables pointing to the various folders of the aforementioned libraries: "$EIGENHOME", "$OPTIMHOME", "$MOSEKHOME" and "$SCSHOME". You can then run:
```bash
make
```
in the root directory of this repo. To run the code, you can then run:
```bash
./run --help
```
to view the help and from there test various solvers applied to the NPA hierarchy.

