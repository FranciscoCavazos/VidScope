// BitcoinNodes.cpp - VidScope Project
// Renders a network of pulsing "Bitcoin nodes" with glowing links and
// an animated transaction dot moving across a path.

#include "../Scenes/Common/CompositeScene.cpp"
#include "../Scenes/Common/CoordinateScene.cpp"
#include "../misc/color.cpp"
#include "../misc/pixels.h"
#include <random>
#include <algorithm>
#include <cmath>

// Bitcoin Nodes Scene - renders the network visualization
class BitcoinNodesScene : public CoordinateScene {
private:
    struct Node {
        float x, y;      // position (normalized 0-1)
        float baseR;     // base radius
        float phase;     // for pulsing
        float hue;       // color hue
        vector<int> nbrs; // neighbor indices
    };
    
    vector<Node> nodes;
    vector<int> transactionPath;
    vector<float> pathSegmentLengths;
    vector<float> pathPrefixSums;
    float totalPathLength;
    mt19937 rng;
    
public:
    BitcoinNodesScene() : CoordinateScene(), rng(42) {
        initializeNodes();
        createConnections();
        createTransactionPath();
    }
    
private:
    void initializeNodes() {
        const int N = 42;
        nodes.reserve(N);
        
        // Layout nodes in 3 clusters
        struct Cluster { float cx, cy, r; } clusters[3] = {
            { 0.30f, 0.45f, 0.12f },  // normalized coordinates
            { 0.60f, 0.40f, 0.14f },
            { 0.50f, 0.70f, 0.13f }
        };
        
        uniform_real_distribution<float> U(0.0f, 1.0f);
        
        for (int i = 0; i < N; ++i) {
            auto &C = clusters[i % 3];
            float ang = U(rng) * 2.0f * M_PI;
            float rad = sqrtf(U(rng)) * C.r;
            float x = C.cx + cosf(ang) * rad;
            float y = C.cy + sinf(ang) * rad;
            
            nodes.push_back({x, y, 0.004f + 0.003f * U(rng), U(rng) * 6.28f, U(rng)});
        }
    }
    
    void createConnections() {
        const int k = 5; // degree
        
        for (int i = 0; i < nodes.size(); ++i) {
            vector<pair<float, int>> distances;
            distances.reserve(nodes.size() - 1);
            
            for (int j = 0; j < nodes.size(); ++j) {
                if (i != j) {
                    float dx = nodes[i].x - nodes[j].x;
                    float dy = nodes[i].y - nodes[j].y;
                    distances.push_back({dx * dx + dy * dy, j});
                }
            }
            
            sort(distances.begin(), distances.end());
            
            for (int t = 0; t < k && t < distances.size(); ++t) {
                nodes[i].nbrs.push_back(distances[t].second);
            }
        }
    }
    
    void createTransactionPath() {
        const int pathLength = 12;
        transactionPath.reserve(pathLength);
        
        int cur = 0;
        vector<int> used(nodes.size(), 0);
        transactionPath.push_back(cur);
        used[cur] = 1;
        
        for (int step = 0; step < pathLength - 1; ++step) {
            auto &nb = nodes[cur].nbrs;
            int best = -1;
            
            // Pick the closest unused neighbor, else any neighbor
            for (int j : nb) {
                if (!used[j]) {
                    best = j;
                    break;
                }
            }
            
            if (best == -1) {
                best = nb[rng() % nb.size()];
            }
            
            cur = best;
            used[cur] = 1;
            transactionPath.push_back(cur);
        }
        
        // Compute path segment lengths
        pathSegmentLengths.reserve(transactionPath.size() - 1);
        pathPrefixSums.push_back(0.0f);
        totalPathLength = 0.0f;
        
        for (size_t i = 0; i + 1 < transactionPath.size(); ++i) {
            auto &a = nodes[transactionPath[i]];
            auto &b = nodes[transactionPath[i + 1]];
            float L = hypotf(a.x - b.x, a.y - b.y);
            pathSegmentLengths.push_back(L);
            totalPathLength += L;
            pathPrefixSums.push_back(totalPathLength);
        }
    }
    
    tuple<float, float, float> hsv2rgb(float h, float s, float v) {
        h = fmodf(h, 1.0f);
        if (h < 0) h += 1.0f;
        
        float i = floorf(h * 6.0f);
        float f = h * 6.0f - i;
        float p = v * (1.0f - s);
        float q = v * (1.0f - f * s);
        float t = v * (1.0f - (1.0f - f) * s);
        
        int ii = int(i) % 6;
        float r, g, b;
        
        if (ii == 0) { r = v; g = t; b = p; }
        else if (ii == 1) { r = q; g = v; b = p; }
        else if (ii == 2) { r = p; g = v; b = t; }
        else if (ii == 3) { r = p; g = q; b = v; }
        else if (ii == 4) { r = t; g = p; b = v; }
        else { r = v; g = p; b = q; }
        
        return {r, g, b};
    }
    
    void drawSoftCircle(float cx, float cy, float rCore, float rGlow,
                       float R, float G, float B, float alphaCore = 0.95f, float alphaGlow = 0.25f) {
        int x0 = max(0, int(floor((cx - rGlow) * VIDEO_WIDTH)));
        int x1 = min(VIDEO_WIDTH - 1, int(ceil((cx + rGlow) * VIDEO_WIDTH)));
        int y0 = max(0, int(floor((cy - rGlow) * VIDEO_HEIGHT)));
        int y1 = min(VIDEO_HEIGHT - 1, int(ceil((cy + rGlow) * VIDEO_HEIGHT)));
        
        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                float px = float(x) / VIDEO_WIDTH;
                float py = float(y) / VIDEO_HEIGHT;
                float dx = px - cx;
                float dy = py - cy;
                float d = sqrtf(dx * dx + dy * dy);
                
                if (d <= rCore) {
                    set_pixel(x, y, R, G, B, alphaCore);
                } else if (d <= rGlow) {
                    float t = (d - rCore) / max(1e-6f, (rGlow - rCore));
                    float a = alphaGlow * (1.0f - t) * (1.0f - t);
                    set_pixel(x, y, R, G, B, a);
                }
            }
        }
    }
    
    void drawGlowLine(float x0, float y0, float x1, float y1,
                      float thickCore, float thickGlow,
                      float R, float G, float B,
                      float alphaCore = 0.5f, float alphaGlow = 0.15f) {
        float dx = x1 - x0;
        float dy = y1 - y0;
        float L = sqrtf(dx * dx + dy * dy);
        if (L < 1e-6f) return;
        
        dx /= L;
        dy /= L;
        
        float pad = thickGlow * 1.2f;
        int minx = max(0, int(floor((min(x0, x1) - pad) * VIDEO_WIDTH)));
        int maxx = min(VIDEO_WIDTH - 1, int(ceil((max(x0, x1) + pad) * VIDEO_WIDTH)));
        int miny = max(0, int(floor((min(y0, y1) - pad) * VIDEO_HEIGHT)));
        int maxy = min(VIDEO_HEIGHT - 1, int(ceil((max(y0, y1) + pad) * VIDEO_HEIGHT)));
        
        for (int y = miny; y <= maxy; ++y) {
            for (int x = minx; x <= maxx; ++x) {
                float px = float(x) / VIDEO_WIDTH;
                float py = float(y) / VIDEO_HEIGHT;
                
                float t = (px - x0) * dx + (py - y0) * dy;
                t = max(0.0f, min(1.0f, t / L));
                
                float cx = x0 + dx * t * L;
                float cy = y0 + dy * t * L;
                float dist = hypotf(px - cx, py - cy);
                
                if (dist <= thickCore) {
                    set_pixel(x, y, R, G, B, alphaCore);
                } else if (dist <= thickGlow) {
                    float u = (dist - thickCore) / max(1e-6f, (thickGlow - thickCore));
                    float a = alphaGlow * (1.0f - u) * (1.0f - u);
                    set_pixel(x, y, R, G, B, a);
                }
            }
        }
    }
    
    void fillVerticalGradient(float r0, float g0, float b0, float r1, float g1, float b1) {
        for (int y = 0; y < VIDEO_HEIGHT; ++y) {
            float t = float(y) / float(VIDEO_HEIGHT - 1);
            float r = r0 * (1 - t) + r1 * t;
            float g = g0 * (1 - t) + g1 * t;
            float b = b0 * (1 - t) + b1 * t;
            
            for (int x = 0; x < VIDEO_WIDTH; ++x) {
                set_pixel(x, y, r, g, b, 1.0f);
            }
        }
    }
    
    void set_pixel(int x, int y, float r, float g, float b, float a) {
        if (x < 0 || y < 0 || x >= VIDEO_WIDTH || y >= VIDEO_HEIGHT) return;
        
        // Simple alpha blending
        int idx = (y * VIDEO_WIDTH + x) * 4;
        float existing_a = pixels[idx + 3];
        float out_a = a + existing_a * (1.0f - a);
        
        if (out_a < 1e-6f) {
            pixels[idx] = pixels[idx + 1] = pixels[idx + 2] = pixels[idx + 3] = 0;
            return;
        }
        
        pixels[idx] = (r * a + pixels[idx] * existing_a * (1.0f - a)) / out_a;
        pixels[idx + 1] = (g * a + pixels[idx + 1] * existing_a * (1.0f - a)) / out_a;
        pixels[idx + 2] = (b * a + pixels[idx + 2] * existing_a * (1.0f - a)) / out_a;
        pixels[idx + 3] = out_a;
    }

public:
    void render() override {
        float t = get_time();
        
        // Clear background with cyber-night gradient
        auto [r0, g0, b0] = hsv2rgb(0.62f, 0.55f, 0.10f);
        auto [r1, g1, b1] = hsv2rgb(0.65f, 0.75f, 0.03f);
        fillVerticalGradient(r0, g0, b0, r1, g1, b1);
        
        // Subtle vignette
        float cx = 0.5f, cy = 0.5f;
        float rad = hypotf(cx, cy) * 1.1f;
        for (int y = 0; y < VIDEO_HEIGHT; ++y) {
            for (int x = 0; x < VIDEO_WIDTH; ++x) {
                float px = float(x) / VIDEO_WIDTH;
                float py = float(y) / VIDEO_HEIGHT;
                float d = hypotf(px - cx, py - cy) / rad;
                float dark = powf(d, 2.2f) * 0.55f;
                set_pixel(x, y, 0, 0, 0, dark);
            }
        }
        
        // Draw links first (glow lines)
        for (int i = 0; i < nodes.size(); ++i) {
            for (int j : nodes[i].nbrs) {
                if (j < i) continue; // avoid double drawing
                
                auto &a = nodes[i];
                auto &b = nodes[j];
                float pul = 0.5f + 0.5f * sinf((a.phase + b.phase) * 0.5f + t * 2.0f);
                float thickCore = 0.0005f + 0.0006f * pul;
                float thickGlow = 0.004f + 0.004f * pul;
                
                auto [r, g, bb] = hsv2rgb(0.55f + 0.05f * pul, 0.6f, 0.9f);
                drawGlowLine(a.x, a.y, b.x, b.y, thickCore, thickGlow, r, g, bb, 0.22f, 0.08f);
            }
        }
        
        // Animated transaction traveling the path
        float travelT = fmodf(t * 0.25f, 1.0f);
        float dist = travelT * totalPathLength;
        
        // Locate segment
        size_t s = 0;
        while (s + 1 < pathPrefixSums.size() && pathPrefixSums[s + 1] < dist) ++s;
        
        if (s < pathSegmentLengths.size()) {
            float segStart = pathPrefixSums[s];
            float u = (dist - segStart) / max(1e-6f, pathSegmentLengths[s]);
            auto &A = nodes[transactionPath[s]];
            auto &B = nodes[transactionPath[s + 1]];
            float x = A.x + (B.x - A.x) * u;
            float y = A.y + (B.y - A.y) * u;
            
            // Transaction glow: bitcoin orange
            auto [r1, g1, b1] = hsv2rgb(0.08f, 0.85f, 1.0f);
            drawSoftCircle(x, y, 0.003f, 0.011f, r1, g1, b1, 0.95f, 0.35f);
            
            // Ring pulses at nodes when dot is near
            if (u < 0.12f || u > 0.88f) {
                auto &N0 = (u < 0.5f) ? A : B;
                float phase = fmodf((u < 0.5f ? u : 1.0f - u) / 0.12f, 1.0f);
                float rad = 0.005f + phase * 0.014f;
                auto [rr, gg, bb] = hsv2rgb(0.08f, 0.8f, 1.0f);
                drawSoftCircle(N0.x, N0.y, rad * 0.65f, rad, rr, gg, bb, 0.35f, 0.18f);
            }
        }
        
        // Draw nodes (pulse + unique hues)
        for (int i = 0; i < nodes.size(); ++i) {
            auto &n = nodes[i];
            float pulse = 0.7f + 0.3f * sinf(n.phase + t * 3.2f + 0.25f * i);
            float rCore = n.baseR * pulse;
            float rGlow = rCore + 0.005f + 0.006f * pulse;
            
            auto [r, g, b] = hsv2rgb(0.45f + 0.15f * n.hue, 0.6f, 1.0f);
            drawSoftCircle(n.x, n.y, rCore, rGlow, r, g, b, 0.95f, 0.28f);
        }
        
        // Subtle scanning sweep
        float sweep = fmodf(t * 0.12f, 1.0f);
        int ySweep = int(sweep * VIDEO_HEIGHT);
        for (int y = max(0, ySweep - 4); y < min(VIDEO_HEIGHT, ySweep + 4); ++y) {
            float a = 0.10f * (1.0f - fabsf(y - ySweep) / 4.0f);
            for (int x = 0; x < VIDEO_WIDTH; ++x) {
                set_pixel(x, y, 0.5f, 0.8f, 1.0f, a);
            }
        }
    }
};

void render_video() {
    // Create the Bitcoin nodes scene
    BitcoinNodesScene bitcoinScene;
    
    // Stage a macroblock with the scene for 12 seconds
    bitcoinScene.stage_macroblock(SilenceBlock(12.0), 1);
    
    // Render the scene
    bitcoinScene.render_microblock();
}
