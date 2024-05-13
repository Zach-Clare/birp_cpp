#include "DataCube.h"

class Camera
{
private:
    int image_dimension;
    std::vector<float> lat;
    std::vector<float> lon;
    int pxn_dist = 200;
    float skysep;
    int ns_lat;
    int ns_lon;
    int nsky;
    std::vector<float> sky1;
    std::vector<float> sky2;
    std::vector<float> xsky;
    std::vector<float> ysky;
    std::vector<std::vector<float>> gei_to_gse;

    void BuildLatLon(float, int);
    void GenerateSky(float);

public:
    DataCube dataCube;
    std::vector<float> image; // image is square, it's size will be image_dimension^2 

    Camera(DataCube cube, float pixel_size_deg, int plot_fov);
    ~Camera();
};