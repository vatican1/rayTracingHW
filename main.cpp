#include <fstream>
#include <iostream>
#include <variant>
#include <vector>
#include <sstream>
#include <optional>

#include <glm/glm.hpp>


inline glm::vec3 rotate(const glm::vec3 v, const glm::vec4 q) {
    glm::vec3 u = glm::vec3(-q.x, -q.y, -q.z);
    float s = q.w;

    return 2.f * dot(u, v) * u +
           (s * s - dot(u, u)) * v +
           2.f * s * cross(u, v);
}

inline glm::vec3 rotate(glm::vec3 v, glm::vec3 axis, float angle)
{
    float sin = std::sin(angle) / 2;
    axis = axis * sin;
    return rotate(v, {axis.x, axis.y, axis.z, std::cos(angle / 2)});
}

struct Ray
{
    Ray(glm::vec3 from_, glm::vec3 direction_) :
        from(from_),
        direction(direction_)
    {}

    glm::vec3 from;
    glm::vec3 direction;
};

struct Plane
{
    glm::vec3 m_normal;
};

struct Ellipsoide
{
    glm::vec3 m_radiuses;
};

struct Box
{
    glm::vec3 m_sizes;
};

struct SceneObject
{
    std::variant<Plane, Ellipsoide, Box> m_primitive;
    glm::vec3 m_position = glm::vec3(0.f, 0.f, 0.f);
    glm::vec4 m_rotation = glm::vec4(0.f, 0.f, 0.f, 1.f);
    glm::vec3 m_color;

    std::optional<float> intersects(Ray ray) const
    {
        ray.from -= m_position;

        ray.from = rotate(ray.from, m_rotation);
        ray.direction = rotate(ray.direction, m_rotation);

        switch (m_primitive.index())
        {
        case 0:
        {
            Plane plane = std::get<Plane>(m_primitive);
            float t = - dot(ray.from, plane.m_normal) / dot(ray.direction, plane.m_normal);

            if (t < 0.f)
                return std::nullopt;
            return std::make_optional(t);
        }
        case 1:
        {
            Ellipsoide ellipsoid = std::get<Ellipsoide>(m_primitive);
            float a, b, c;
            glm::vec3 o_div_r = ray.from      / ellipsoid.m_radiuses;
            glm::vec3 d_div_r = ray.direction / ellipsoid.m_radiuses;
            a = dot(d_div_r, d_div_r);
            b = dot(o_div_r, d_div_r);
            c = dot(o_div_r, o_div_r);

            float h = b * b - a * (c - 1.f);
            if( h < 0.f )
                return std::nullopt;

            h = std::sqrt(h);
            float t1 = (-b - h) / a;
            float t2 = (-b + h) / a;

            if (t2 < 0)
                return std::nullopt;

            if (t1 < 0)
                t1 = t2;
            return std::make_optional(t1);
        }
        case 2:
        {
            Box box = std::get<Box>(m_primitive);
            glm::vec3 t1;
            glm::vec3 t2;
            for (int i = 0; i < 3; ++i)
            {
                t1[i] = (-box.m_sizes[i] - ray.from[i]) / ray.direction[i];
                t2[i] = (box.m_sizes[i] - ray.from[i]) / ray.direction[i];
                if (t1[i] > t2[i])
                    std::swap(t1[i], t2[i]);
            }

            float t1_max = std::max(std::max(t1.x, t1.y), t1.z);
            float t2_min = std::min(std::min(t2.x, t2.y), t2.z);

            if (t1_max > t2_min || t2_min < 0)
                return std::nullopt;

            if (t1_max < 0)
                t1_max = t2_min;
            return std::make_optional(t1_max);
        }
        default:
            break;
        }
        return std::nullopt;
    }
};

struct OutputColor
{
    char r;
    char g;
    char b;
};

struct Canvas
{

    Canvas(size_t w_, size_t h_):
        w(w_),
        h(h_),
        data(w_ * h_)
    {}

    std::vector<OutputColor> data;
    size_t w;
    size_t h;

    void setPixelColor(const size_t i, const size_t j, const glm::vec3 & colorF)
    {
        OutputColor color;
        color.r = colorF.x * 255.0;
        color.g = colorF.y * 255.0;
        color.b = colorF.z * 255.0;
        data.at(i + j * w) = color;
    }

    void save(const std::string & fileName) const
    {
        std::ofstream out(fileName, std::ios_base::binary);

        if (!out)
            throw std::runtime_error("File open error");

        out << "P6\n"
            << w << " " << h << std::endl
            << "255" << std::endl;

        char *data_char = (char *)(void *)data.data();
        out.write(data_char, data.size() * sizeof(OutputColor));

        out.close();
    }
};

struct CameraParams
{
    glm::vec3 m_cameraPosition;

    glm::vec3 m_cameraRight;
    glm::vec3 m_cameraUp;
    glm::vec3 m_cameraForward;

    float m_cameraFovX;
    float m_cameraFovY = 0.f;

    void setFovY(int w, int h)
    {
        m_cameraFovY = 2.f * std::atan(std::tan(m_cameraFovX * 0.5f) * float(h) / float(w));
    }
};

struct Camera
{
public:

    Camera(const CameraParams & params):
        m_params(params)
    {}

    glm::vec3 getXYZVectors(int i)
    {
        switch (i) {
        case 0:
            return m_params.m_cameraRight;
        case 1:
            return m_params.m_cameraUp;
        case 2:
            return m_params.m_cameraForward;
        default:
            return glm::vec3(0.f, 0.f, 0.f);
        }
    }

    glm::vec3 getPosition()
    {
        return m_params.m_cameraPosition;
    }

    float getFovY()
    {
        assert(m_params.m_cameraFovY != 0.f);
        return m_params.m_cameraFovY;
    }
    float getFovX()
    {
        return m_params.m_cameraFovX;
    }

private:
    CameraParams m_params;
};


struct Scene
{
    glm::vec3 m_bgColor;

    Canvas * p_canvas;
    Camera * p_camera;

    std::vector<SceneObject> m_sceneObjects;

    Ray canvasToCamera(int x, int y) const
    {
        glm::vec3 t;
        t.x =  (2.f * (float(x) + 0.5f) / float(p_canvas->w) - 1.0) * std::tan(p_camera->getFovX() / 2);
        t.y = -(2.f * (float(y) + 0.5f) / float(p_canvas->h) - 1.0) * std::tan(p_camera->getFovY() / 2);
        t.z = 1;

        glm::vec3 d(0.f, 0.f, 0.f);
        for (int i = 0; i < 3; ++i) {
            d += t[i] * p_camera->getXYZVectors(i);
        }

        return Ray(p_camera->getPosition(), d);
    }

    void render() const
    {
        if(!p_canvas)
        {
            std::cout << "canvas not created\n";
            return;
        }

        for (int i = 0; i < p_canvas->w; ++i)
        {
            for (int j = 0; j < p_canvas->h; ++j)
            {
                Ray ray = canvasToCamera(i, j);

                glm::vec3 color = m_bgColor;
                float distance = 1e9;

                for (const SceneObject & object : m_sceneObjects)
                {
                    std::optional<float> intersection = object.intersects(ray);
                    if (intersection && intersection.value() < distance)
                    {
                        distance = intersection.value();
                        color = object.m_color;
                    }
                }

                p_canvas->setPixelColor(i, j, color);
            }
        }
    }


    void save(const std::string & fileName) const
    {
        p_canvas->save(fileName);
    }
};


const Scene * parceScene(const std::string & fileName)
{
    std::ifstream file(fileName);
    if (!file)
    {
        std::cout << "ifstream not opened\n";
        return nullptr;
    }

    Scene * p_scene = new Scene;
    CameraParams cameraParams;
    size_t canvasW, canvasH;

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "DIMENSIONS") {
            iss >> canvasW >> canvasH;
        } else if (key == "BG_COLOR") {
            iss >> p_scene->m_bgColor.r >> p_scene->m_bgColor.g >> p_scene->m_bgColor.b;
        } else if (key == "CAMERA_POSITION") {
            iss >> cameraParams.m_cameraPosition.x >> cameraParams.m_cameraPosition.y >> cameraParams.m_cameraPosition.z;
        } else if (key == "CAMERA_RIGHT") {
            iss >> cameraParams.m_cameraRight.x >> cameraParams.m_cameraRight.y >> cameraParams.m_cameraRight.z;
        } else if (key == "CAMERA_UP") {
            iss >> cameraParams.m_cameraUp.x >> cameraParams.m_cameraUp.y >> cameraParams.m_cameraUp.z;
        } else if (key == "CAMERA_FORWARD") {
            iss >> cameraParams.m_cameraForward.x >> cameraParams.m_cameraForward.y >> cameraParams.m_cameraForward.z;
        } else if (key == "CAMERA_FOV_X") {
            iss >> cameraParams.m_cameraFovX;
        } else if(key == "NEW_PRIMITIVE") {
            break;
        } else {
            continue;
        }
    }
    p_scene->p_canvas = new Canvas(canvasW, canvasH);
    cameraParams.setFovY(canvasW, canvasH);
    p_scene->p_camera = new Camera(cameraParams);

    SceneObject object;
    bool isNewPrimitive = false;
    do
    {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "NEW_PRIMITIVE")
        {
            if (isNewPrimitive) {
                p_scene->m_sceneObjects.push_back(object);
            }
            object = SceneObject();
            isNewPrimitive = true;
        } else if (key == "ELLIPSOID" || key == "PLANE" || key == "BOX") {
            if(key == "ELLIPSOID")
            {
                Ellipsoide el;
                iss >> el.m_radiuses.x >> el.m_radiuses.y >> el.m_radiuses.z;
                object.m_primitive = el;
            }
            else if(key == "BOX")
            {
                Box box;
                iss >> box.m_sizes.x >> box.m_sizes.y >> box.m_sizes.z;
                object.m_primitive = box;
            }
            else if(key == "PLANE")
            {
                Plane plane;
                iss >> plane.m_normal.x >> plane.m_normal.y >> plane.m_normal.z;
                object.m_primitive = plane;
            }

        } else if (key == "POSITION") {
            iss >> object.m_position.x >> object.m_position.y >> object.m_position.z;
        } else if (key == "COLOR") {
            iss >> object.m_color.x >> object.m_color.y >> object.m_color.z;
        } else if (key == "ROTATION") {
            iss >> object.m_rotation.x >> object.m_rotation.y
                >> object.m_rotation.z >> object.m_rotation.w;
        }
    } while (std::getline(file, line));
    p_scene->m_sceneObjects.push_back(object);
    file.close();

    return p_scene;
}


int main(int argc, char *argv[])
{
    assert(argc == 3);
    const Scene * scene = parceScene(argv[1]);
    scene->render();
    scene->save(argv[2]);
    return 0;
}

