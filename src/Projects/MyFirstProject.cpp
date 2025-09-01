#include "../Scenes/Media/LatexScene.cpp"
#include "../Scenes/Common/CompositeScene.cpp"

void render_video() {
    // Create a simple LaTeX scene
    LatexScene latex("Hello SwapTube!", 0.5);
    
    // Stage a macroblock (audio segment) with 1 microblock (video segment)
    latex.stage_macroblock(SilenceBlock(3.0), 1);
    
    // Render the microblock
    latex.render_microblock();
    
    // Add another scene with a transition
    latex.begin_latex_transition(MICRO, "Welcome to SwapTube!");
    latex.stage_macroblock(SilenceBlock(2.0), 1);
    latex.render_microblock();
}