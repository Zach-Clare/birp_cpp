#include "DataCube.h"

class Camera
{
private:
    int image_dimension;
    std::vector<float> lat;
    std::vector<float> lon;
    int fov;
    int pxn_dist = 200;
    float skysep;
    int ns_lat;
    int ns_lon;
    int nsky;
    std::vector<std::vector<float>> gei_to_gse;
    struct {
        float x;
        float y;
        float z;
    } position;
    struct {
        float x;
        float y;
        float z;
    } aim;
    std::vector<float> ray_dist;
    std::vector<float> ray_widths;
    float pixel_size_deg;
    float pixel_size_rad;
    float ray_width;
    int ray_samples;

    void BuildLatLon(float, int);
    void GenerateSky(float);
    void GetSky();
    void GenerateRayDistWidth(int);
    void Integrate();

    float Orient();
    std::vector<float> PointToPlane(std::vector<float>, float, float, std::vector<float> distance_vec, std::vector<float> unit_vec, std::vector<float> north, std::vector<float> right);

public:
    DataCube dataCube;
    float image[144][144] {0}; // image is square, it's size will be image_dimension^2 
    
    std::vector<float> sky1;
    std::vector<float> sky2;
    std::vector<float> xsky;
    std::vector<float> ysky;

    Camera(DataCube cube, float pixel_size_deg, int plot_fov);
    ~Camera();

    void SetPosition(float, float, float);
    void SetAim(float, float, float);
    void Render();
    int ToDat(std::string filename);
    int ToFITS(std::string filename);
};