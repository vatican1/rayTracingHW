#include <fstream>
#include <iostream>
#include <variant>
#include <vector>
#include <sstream>
#include <optional>

#include <glm/glm.hpp>

#include <random>


inline glm::vec3 rotate(const glm::vec3 v, const glm::vec4 q) {
    const glm::vec3 u = glm::vec3(q.x, q.y, q.z);
    float s = q.w;

    return 2.f * dot(u, v) * u +
           (s * s - dot(u, u)) * v +
           2.f * s * cross(u, v);
}

inline glm::vec3 saturate(const glm::vec3 &color) {
    return clamp(color, {0.f, 0.f, 0.f}, {1.f, 1.f, 1.f});
}

inline glm::vec3 acesTonemap(const glm::vec3 & x) {
    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;
    return saturate((x*(a*x+b))/(x*(c*x+d)+e));
}


glm::vec3 generateSphereDir()
{

    std::uniform_real_distribution<> genH(-1.f, 1.f);
    std::uniform_real_distribution<> genAng(0, 2.f * M_PI);
    static std::minstd_rand gen;

    const float ang = genAng(gen);
    const float height = genH(gen);

    const float proj = std::sqrt(1 - height * height);

    return glm::vec3(proj * std::cos(ang), proj * std::sin(ang), height);
}

struct Intersection
{
    float distance;
    glm::vec3 normal;
    bool fromOutside;

    glm::vec3 color{0.f, 0.f, 0.f};
};

struct Ray
{
    Ray(glm::vec3 from_, glm::vec3 direction_) :
        from(from_),
        direction(glm::normalize(direction_))
    {}

    glm::vec3 from;
    glm::vec3 direction;

    void shiftStart(const float & dist = 2e-3)
    {
        from += direction * dist;
    }
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


float t_Plane(const Plane & plane, const Ray & ray)
{
    return - dot(ray.from, plane.m_normal) / dot(ray.direction, plane.m_normal);
}

std::pair<float, float> t1_t2_Ellipsiod(const Ellipsoide & ellipsoid, const Ray & ray, bool * ok)
{
    *ok = true;
    float a, b, c;
    glm::vec3 o_div_r = ray.from      / ellipsoid.m_radiuses;
    glm::vec3 d_div_r = ray.direction / ellipsoid.m_radiuses;
    a = dot(d_div_r, d_div_r);
    b = dot(o_div_r, d_div_r);
    c = dot(o_div_r, o_div_r);

    float h = b * b - a * (c - 1.f);
    if( h < 0.f )
    {
        *ok = false;
        std::make_pair(0.f, 0.f);
    }
    h = std::sqrt(h);
    float t1 = (-b - h) / a;
    float t2 = (-b + h) / a;

    return std::make_pair(t1, t2);
}

std::pair<float, float> t1Max_t2Min_Box(const Box & box, const Ray & ray)
{
    glm::vec3 t1, t2;
    for (int i = 0; i < 3; ++i)
    {
        t1[i] = (-box.m_sizes[i] - ray.from[i]) / ray.direction[i];
        t2[i] = (box.m_sizes[i] - ray.from[i]) / ray.direction[i];
        if (t1[i] > t2[i])
            std::swap(t1[i], t2[i]);
    }

    float t1Max = std::max(std::max(t1.x, t1.y), t1.z);
    float t2Min = std::min(std::min(t2.x, t2.y), t2.z);

    return std::make_pair(t1Max, t2Min);
}

enum OBJECT_MATERIAL
{
    OM_METALLIC = 0,
    OM_DIELECTRIC,
    OM_DIFFUSE
};


struct SceneObject
{
    std::variant<Plane, Ellipsoide, Box> m_primitive;
    glm::vec3 m_position = glm::vec3(0.f, 0.f, 0.f);
    glm::vec4 m_rotation = glm::vec4(0.f, 0.f, 0.f, 1.f);
    glm::vec3 m_color;
    glm::vec3 m_emission = glm::vec3(0.f, 0.f, 0.f);;

    OBJECT_MATERIAL m_objectMaterial = OM_DIFFUSE;
    float m_ior;


    void prepareRay(Ray & ray) const
    {
        ray.from -= m_position;

        ray.from = rotate(ray.from, {-m_rotation.x,-m_rotation.y,-m_rotation.z, m_rotation.w} );
        ray.direction = rotate(ray.direction, {-m_rotation.x,-m_rotation.y,-m_rotation.z, m_rotation.w});
    }

    // Вызывается, если уверенны, что пересечение есть!
    Intersection intersectFull(Ray ray) const
    {
        prepareRay(ray);
        switch (m_primitive.index())
        {
        case 0:
        {
            const Plane & plane = std::get<Plane>(m_primitive);
            float t = t_Plane (plane, ray);
            assert(t >= 0.f);

            Intersection intersection;
            intersection.distance = t;
            intersection.normal = plane.m_normal;
            intersection.fromOutside = true;

            intersection.normal = rotate(intersection.normal, m_rotation);

            return intersection;
        }
        case 1:
        {
            Ellipsoide ellipsoid = std::get<Ellipsoide>(m_primitive);

            bool ok;
            float t1, t2;
            std::tie(t1, t2) = t1_t2_Ellipsiod(ellipsoid, ray, &ok);
            assert(ok);

            Intersection intersection;
            intersection.distance = t1 > 0 ? t1 : t2;
            intersection.fromOutside = t1 > 0;
            glm::vec3 interPoint = ray.from + intersection.distance * ray.direction;
            intersection.normal = glm::normalize( interPoint / ellipsoid.m_radiuses / ellipsoid.m_radiuses);
            intersection.normal = rotate(intersection.normal, m_rotation);

            if (!intersection.fromOutside)
                intersection.normal = -intersection.normal;

            return intersection;

        }
        case 2:
        {
            Box box = std::get<Box>(m_primitive);

            float t1Max, t2Min;

            std::tie(t1Max, t2Min) = t1Max_t2Min_Box(box, ray);

            // оставлено на дальнейшее растаскивание
            Intersection intersection;
            intersection.distance = t1Max > 0 ? t1Max : t2Min;
            intersection.fromOutside = t1Max > 0;
            glm::vec3 interPoint = ray.from + intersection.distance * ray.direction;

            intersection.normal = interPoint / box.m_sizes;

            int maxIdx = 0;
            {
                float maxDist = 0.f;
                for (int i = 0; i < 3; ++i)
                {
                    if (std::abs(intersection.normal[i]) >= maxDist)
                    {
                        maxDist = std::abs(intersection.normal[i]);
                        maxIdx = i;
                    }
                }
            }

            intersection.normal = {(maxIdx == 0 ? (intersection.normal[0] > 0.f ? 1.f : -1.f) : 0.f),
                                   (maxIdx == 1 ? (intersection.normal[1] > 0.f ? 1.f : -1.f) : 0.0),
                                   (maxIdx == 2 ? (intersection.normal[2] > 0.f ? 1.f : -1.f) : 0.0)};

            if (!intersection.fromOutside)
                intersection.normal = -intersection.normal;

            intersection.normal = rotate(intersection.normal, m_rotation);

            return intersection;
        }
        default:
            break;
        }
        return Intersection();
    }

    // Возвращает либо больше нуля, либо std::nullopt
    std::optional<float> intersectsSimple(Ray ray) const
    {
        prepareRay(ray);
        switch (m_primitive.index())
        {
        case 0:
        {
            const Plane & plane = std::get<Plane>(m_primitive);
            float t = t_Plane (plane, ray);

            if (t < 0.f)
                return std::nullopt;
            return std::make_optional(t);
        }
        case 1:
        {
            Ellipsoide ellipsoid = std::get<Ellipsoide>(m_primitive);
            bool ok;
            float t1, t2;
            std::tie(t1, t2) = t1_t2_Ellipsiod(ellipsoid, ray, &ok);
            if(!ok)
                return std::nullopt;

            if (t2 < 0)
                return std::nullopt;

            if (t1 < 0)
                t1 = t2;
            return std::make_optional(t1);
        }
        case 2:
        {
            Box box = std::get<Box>(m_primitive);

            float t1Max, t2Min;

            std::tie(t1Max, t2Min) = t1Max_t2Min_Box(box, ray);

            if (t1Max > t2Min || t2Min < 0)
                return std::nullopt;

            if (t1Max < 0)
                t1Max = t2Min;
            return std::make_optional(t1Max);
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

    inline void setPixelColor(const size_t i, const size_t j, const glm::vec3 & colorF)
    {
        OutputColor color;
        color.r = std::round(std::clamp(colorF.x * 255.0, 0.0, 255.0));
        color.g = std::round(std::clamp(colorF.y * 255.0, 0.0, 255.0));
        color.b = std::round(std::clamp(colorF.z * 255.0, 0.0, 255.0));
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

struct Light
{
    explicit Light(const glm::vec3 & m_lightIntensity)
        : m_lightIntensity(m_lightIntensity)
    {}

    virtual Ray rayFrom(glm::vec3 from) const = 0;
    virtual bool isObscuredObject(const Ray & ray, const float distToObj) const = 0;
    //свет из источника в точке
    virtual glm::vec3 colorDist(const glm::vec3 & from) const = 0;

    glm::vec3 m_lightIntensity;
};

struct PointLight : public Light
{
    PointLight(const glm::vec3 & lightIntensity,
               const glm::vec3 & lightPosition,
               const glm::vec3 & lightAttenuation):
        m_lightPosition(lightPosition),
        m_lightAttenuation(lightAttenuation),
        ::Light(lightIntensity)
    {}

    virtual Ray rayFrom(glm::vec3 from) const
    {
        Ray ray(from, m_lightPosition - from);
        return ray;
    }

    virtual bool isObscuredObject(const Ray & ray, const float distToObj) const
    {
        return glm::length(ray.from - m_lightPosition) > distToObj;
    }


    virtual glm::vec3 colorDist(const glm::vec3 & from) const
    {
        float dist = glm::length(from - m_lightPosition);
        return m_lightIntensity * (1.f / (m_lightAttenuation[0] +
                                          m_lightAttenuation[1] * dist +
                                          m_lightAttenuation[2] * dist * dist));
    }


    glm::vec3 m_lightPosition;
    glm::vec3 m_lightAttenuation;
};

struct DirectionalLight : public Light
{
    DirectionalLight(const glm::vec3 & lightIntensity,
                     const glm::vec3 & lightDirection):
        m_lightDirection(lightDirection),
        ::Light(lightIntensity)
    {}

    virtual Ray rayFrom(glm::vec3 from) const
    {
        Ray ray(from, m_lightDirection);
        return ray;
    }

    virtual bool isObscuredObject(const Ray & ray, const float distToObj) const
    {
        return distToObj > 0.f;
    }

    virtual glm::vec3 colorDist(const glm::vec3 & from) const
    {
        return m_lightIntensity;
    }

    glm::vec3 m_lightDirection;
};

enum LIGHT_TYPE
{
    LT_UNKNOWN,
    LT_POINT,
    LT_DIRECTIONAL
};


struct Scene
{
    glm::vec3 m_bgColor;

    int m_rayDepth;
    int m_samples;
    glm::vec3 m_ambientLight;

    Canvas * p_canvas;
    Camera * p_camera;

    std::vector<Light *> m_lights;

    std::vector<SceneObject> m_sceneObjects;

    Ray canvasToCamera(int x, int y) const
    {
        std::uniform_real_distribution<> genShift(0.f, 1.f);
        static std::minstd_rand gen;
        float shift = genShift(gen);

        glm::vec3 t;
        t.x =  (2.f * (float(x) + shift) / float(p_canvas->w) - 1.0) * std::tan(p_camera->getFovX() / 2);
        t.y = -(2.f * (float(y) + shift) / float(p_canvas->h) - 1.0) * std::tan(p_camera->getFovY() / 2);
        t.z = 1;

        glm::vec3 d(0.f, 0.f, 0.f);
        for (int i = 0; i < 3; ++i) {
            d += t[i] * p_camera->getXYZVectors(i);
        }

        return Ray(p_camera->getPosition(), d);
    }

    glm::vec3 colorFromRay(const Ray & ray, int currDepth) const
    {
        ++currDepth;
        if(currDepth > m_rayDepth)
        {
            return glm::vec3(0.f, 0.f, 0.f);
        }

        int objNumber;
        float distanceToObject;
        std::tie(objNumber, distanceToObject) = findNearestObject(ray);
        if(objNumber == -1)
        {
            return m_bgColor;
        }

        const SceneObject & object = m_sceneObjects[objNumber];
        Intersection intersection = object.intersectFull(ray);
        intersection.color = object.m_emission;

        glm::vec3 intersectionEnd = ray.from + ray.direction * intersection.distance; // Пересечение объекта и луча из камеры

        if(object.m_objectMaterial == OM_DIFFUSE)
        {
            glm::vec3 direction = generateSphereDir();
            if (dot(direction, intersection.normal) < 0)
                direction = -direction;
            Ray reflectedRay (intersectionEnd, direction);
            reflectedRay.shiftStart();
            const glm::vec3 reflectedColor = colorFromRay(reflectedRay, currDepth);

            intersection.color += reflectedColor * object.m_color * 2.f * dot(direction, intersection.normal);
            return  intersection.color;
        }
        else if(object.m_objectMaterial == OM_DIELECTRIC)
        {
            const float cosIn = -dot(intersection.normal, ray.direction);
            const float eta1 = intersection.fromOutside ? 1.f : object.m_ior;
            const float eta2 = intersection.fromOutside ? object.m_ior : 1.f;

            const float sinOut = (eta1 / eta2) * std::sqrt(1 - cosIn * cosIn);

            float reflectionCoefficient = 1.f;
            float randomDir = 0.f;
            if (sinOut < 1.f)
            {
                const float R0 = std::pow((eta1 - eta2) / (eta1 + eta2), 2.f);
                reflectionCoefficient = R0 + (1.f - R0) * std::pow(1.f - cosIn, 5.f);
                std::uniform_real_distribution<> genDir(0.f, 1.f);
                static std::minstd_rand gen;
                randomDir = genDir(gen);
            }
//!---------------------------------------------------------------------------------------------------
            if(randomDir < reflectionCoefficient)
            {
                glm::vec3 dir1 = ray.direction - intersection.normal * dot(intersection.normal, ray.direction) * 2.f;
                Ray reflectRay(intersectionEnd, dir1);
                reflectRay.shiftStart();
                const glm::vec3 appendColor1 = colorFromRay(reflectRay, currDepth);

                intersection.color += appendColor1;

            }
//!---------------------------------------------------------------------------------------------------
            else
            {
                if(sinOut >= 1.f)
                {
                    assert(false);
                    return intersection.color;
                }

                float cosOut = std::sqrt(1 - sinOut * sinOut);
                const glm::vec3 dir2 = (eta1 / eta2) * ray.direction + ((eta1 / eta2) * cosIn - cosOut) * intersection.normal;

                Ray reflectRay2(intersectionEnd, dir2);
                reflectRay2.shiftStart();
                const glm::vec3 appendColor2 = colorFromRay(reflectRay2, currDepth);

                glm::vec3 currColor(0.f, 0.f, 0.f);

                currColor = appendColor2;
                if(intersection.fromOutside)
                    currColor *= object.m_color;

                intersection.color += currColor;
            }
//!---------------------------------------------------------------------------------------------------
            return intersection.color;
        }
        else if(object.m_objectMaterial == OM_METALLIC)
        {
            glm::vec3 dir = ray.direction - intersection.normal * dot(intersection.normal, ray.direction) * 2.f;
            Ray reflectRay(intersectionEnd, dir);
            reflectRay.shiftStart();
            const glm::vec3 appendColor = colorFromRay(reflectRay, currDepth);

            intersection.color += object.m_color * appendColor;


            return intersection.color;
        }
        assert(false);
        return glm::vec3(0.f, 0.f, 0.f);
    }

    void render() const
    {
        if(!p_canvas)
        {
            std::cout << "canvas not created\n";
            return;
        }

#pragma  omp parallel for collapse(2)
        for (int i = 0; i < p_canvas->w; ++i)
        {
            for (int j = 0; j < p_canvas->h; ++j)
            {
                glm::vec3 color(0.f, 0.f, 0.f);
                for(int s = 0; s < m_samples; ++s)
                {
                    const Ray ray = canvasToCamera(i, j);
                    color += colorFromRay(ray, 0);
                }
                color /= float(m_samples);
                color = acesTonemap(color);
                color = glm::pow(color, glm::vec3(1.f / 2.2f));

                p_canvas->setPixelColor(i, j, color);

            }
        }
    }

    // Возвращает номер ближайшего по лучу объекта или -1 и дистанцию
    std::pair<int, float> findNearestObject(const Ray & ray) const
    {
        float distance = 1e9;
        int objectNumber = -1;
        for (size_t i = 0; i < m_sceneObjects.size(); ++i)
        {
            const SceneObject & object = m_sceneObjects[i];
            std::optional<float> intersection = object.intersectsSimple(ray);
            if (intersection && intersection.value() < distance)
            {
                distance = intersection.value();
                objectNumber = int(i);
            }
        }
        return std::make_pair(objectNumber, objectNumber == -1 ? -1.f : distance);
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

    LIGHT_TYPE lightType = LT_UNKNOWN;

    // Все возможные параметры источников
    glm::vec3 lightIntensity;

    glm::vec3 lightPosition;
    glm::vec3 lightAttenuation;

    glm::vec3 lightDirection;

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
            float x, y, z;
            iss >> x >> y >> z;
            cameraParams.m_cameraPosition = glm::vec3(x, y, z);
        } else if (key == "CAMERA_RIGHT") {
            float x, y, z;
            iss >> x >> y >> z;
            cameraParams.m_cameraRight = glm::vec3(x, y, z);
        } else if (key == "CAMERA_UP") {
            float x, y, z;
            iss >> x >> y >> z;
            cameraParams.m_cameraUp = glm::vec3(x, y ,z);
        } else if (key == "CAMERA_FORWARD") {
            float x, y, z;
            iss >> x >> y >> z;
            cameraParams.m_cameraForward = glm::vec3(x, y, z);
        } else if (key == "CAMERA_FOV_X") {
            iss >> cameraParams.m_cameraFovX;
        } else if (key == "RAY_DEPTH") {
            iss >> p_scene->m_rayDepth;
        } else if (key == "SAMPLES") {
            iss >> p_scene->m_samples;
        } else if (key == "AMBIENT_LIGHT") {
            float x, y, z;
            iss >> x >> y >> z;
            p_scene->m_ambientLight = glm::vec3(x, y, z);
        } else if(key == "NEW_LIGHT") {
            if(lightType != LT_UNKNOWN)
            {
                if(lightType == LT_POINT)
                {
                    p_scene->m_lights.push_back(new PointLight(lightIntensity, lightPosition, lightAttenuation));
                }
                else if( lightType == LT_DIRECTIONAL)
                {
                    p_scene->m_lights.push_back(new DirectionalLight(lightIntensity, lightDirection));
                }
            }
        } else if (key == "LIGHT_INTENSITY") {
            float x, y, z;
            iss >> x >> y >> z;
            lightIntensity = {x, y ,z};
        } else if (key == "LIGHT_DIRECTION") {
            lightType = LT_DIRECTIONAL;
            float x, y, z;
            iss >> x >> y >> z;
            lightDirection = {x, y ,z};
        } else if (key == "LIGHT_POSITION") {
            lightType = LT_POINT;
            float x, y, z;
            iss >> x >> y >> z;
            lightPosition = {x, y ,z};
        } else if (key == "LIGHT_ATTENUATION") {
            lightType = LT_POINT;
            float x, y, z;
            iss >> x >> y >> z;
            lightAttenuation = {x, y, z};
        } else if(key == "NEW_PRIMITIVE") {
            break;
        } else
        {
            continue;
        }
    }


    if(lightType == LT_POINT)
    {
        p_scene->m_lights.push_back(new PointLight(lightIntensity, lightPosition, lightAttenuation));
    }
    else if( lightType == LT_DIRECTIONAL)
    {
        p_scene->m_lights.push_back(new DirectionalLight(lightIntensity, lightDirection));
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
        } else if (key == "EMISSION") {
            iss >> object.m_emission.x >> object.m_emission.y >> object.m_emission.z;
        } else if (key == "ROTATION") {
            iss >> object.m_rotation.x >> object.m_rotation.y
                >> object.m_rotation.z >> object.m_rotation.w;
        }
        else if (key == "METALLIC" || key == "DIELECTRIC") {
            if(key == "METALLIC" )
                object.m_objectMaterial = OM_METALLIC;
            else if(key == "DIELECTRIC" )
                object.m_objectMaterial = OM_DIELECTRIC;
        }
        else if (key == "IOR") {
            iss >> object.m_ior;
        }
    } while (std::getline(file, line));
    p_scene->m_sceneObjects.push_back(object);
    file.close();

    return p_scene;
}


int main(int argc, char *argv[])
{
    std::cout << "start" << std::endl;
    assert(argc == 3);
    const Scene * scene = parceScene(argv[1]);
    std::cout << "end parce" << std::endl;
    scene->render();
    std::cout << "end render" << std::endl;
    scene->save(argv[2]);
    std::cout << "finish" << std::endl;
    return 0;
}

