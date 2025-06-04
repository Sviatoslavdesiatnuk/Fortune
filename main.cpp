#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <memory>
#include "voronoi.h"

namespace {
    std::unique_ptr<Voronoi> voronoiEngine;
    std::vector<std::unique_ptr<VoronoiPoint>> sites;
    std::vector<VEdge> edges;

    void initGL() {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glPointSize(5.0);
        glLineWidth(1.0);
    }

    void display() {
        glClear(GL_COLOR_BUFFER_BIT);

        glColor3f(0.0f, 1.0f, 1.0f);
        glBegin(GL_POINTS);
        for (const auto& site : sites)
            glVertex2d(site->x, site->y);
        glEnd();

        glColor3f(0.0f, 0.8f, 0.5f);
        glBegin(GL_LINES);
        for (const auto& edge : edges) {
            glVertex2d(edge.VertexA.x, edge.VertexA.y);
            glVertex2d(edge.VertexB.x, edge.VertexB.y);
        }
        glEnd();

        glFlush();
    }

    void reshape(GLsizei width, GLsizei height) {
        if (height == 0) height = 1;
        const GLfloat aspect = static_cast<GLfloat>(width) / height;

        glViewport(0, 0, width, height);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        if (width >= height)
            gluOrtho2D(-3.0 * aspect, 3.0 * aspect, -3.0, 3.0);
        else
            gluOrtho2D(-3.0, 3.0, -3.0 / aspect, 3.0 / aspect);
    }

    void prepareSites() {
        const std::vector<std::pair<double, double>> coords = {
            {0.40, 0.730}, {0.880, 0.870}, {0.550, 0.430}, {0.560, 0.750},
            {0.460, 0.180}, {0.70, 0.110}, {0.380, 0.920}, {0.660, 0.30},
            {0.800, 0.470}, {0.190, 0.580}, {0.960, 0.30}, {0.270, 0.350}
        };

        for (const auto& [x, y] : coords)
            sites.push_back(std::make_unique<VoronoiPoint>(x, y));
    }

    void runVoronoi() {
        voronoiEngine = std::make_unique<Voronoi>();

        std::vector<VoronoiPoint*> rawSites;
        rawSites.reserve(sites.size());
        for (const auto& site : sites)
            rawSites.push_back(site.get());

        constexpr double minY = -10.0;
        constexpr double maxY = 10.0;
        edges = voronoiEngine->ComputeVoronoiGraph(rawSites, minY, maxY);

        std::cout << "\t\tVoronoi Edges\n";
        for (const auto& edge : edges) {
            std::cout << "(" << edge.VertexA.x << ", " << edge.VertexA.y << ")\t(";
            std::cout << edge.VertexB.x << ", " << edge.VertexB.y << ")\n";
        }

        for (const auto& site : sites) {
            std::cout << "\t\tSite =  (" << site->x << ", " << site->y << ")\n";
            for (const auto& edge : edges) {
                if ((edge.Left_Site.x == site->x && edge.Left_Site.y == site->y) ||
                    (edge.Right_Site.x == site->x && edge.Right_Site.y == site->y)) {
                    std::cout << "(" << edge.VertexA.x << ", " << edge.VertexA.y << ")\t(";
                    std::cout << edge.VertexB.x << ", " << edge.VertexB.y << ")\n";
                }
            }
            std::cout << "\n";
        }
    }
}

int main(int argc, char** argv) {
    system("mode con: cols=100 lines=2500");

    prepareSites();
    runVoronoi();

    glutInit(&argc, argv);
    glutInitWindowSize(840, 768);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Voronoi Diagram");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    initGL();
    glutMainLoop();

    return 0;
}