/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifdef _WIN32

#include "wscene/canvas/win_opengl_canvas.h"
#include "GL/glew.h"

namespace wgp {

	WinOpenGLCanvas::WinOpenGLCanvas(HWND hWnd) {
        m_hWnd = hWnd;
        m_hrc = NULL;
        HDC hdc = GetDC(hWnd);
        if (hdc) {
            PIXELFORMATDESCRIPTOR pfd =											// pfd Tells Windows How We Want Things To Be
            {
                sizeof(PIXELFORMATDESCRIPTOR),									// Size Of This Pixel Format Descriptor
                1,																// Version Number
                PFD_DRAW_TO_WINDOW |											// Format Must Support Window
                PFD_SUPPORT_OPENGL |											// Format Must Support OpenGL
                PFD_DOUBLEBUFFER,												// Must Support Double Buffering
                PFD_TYPE_RGBA,													// Request An RGBA Format
                24,										                        // Select Our Color Depth
                0, 0, 0, 0, 0, 0,												// Color Bits Ignored
                0,																// No Alpha Buffer
                0,																// Shift Bit Ignored
                0,																// No Accumulation Buffer
                0, 0, 0, 0,														// Accumulation Bits Ignored
                16,																// 16Bit Z-Buffer (Depth Buffer)  
                0,																// No Stencil Buffer
                0,																// No Auxiliary Buffer
                PFD_MAIN_PLANE,													// Main Drawing Layer
                0,																// Reserved
                0, 0, 0															// Layer Masks Ignored
            };
            int pixel_format = ChoosePixelFormat(hdc, &pfd);
            if (pixel_format) {
                if (SetPixelFormat(hdc, pixel_format, &pfd)) {
                    m_hrc = wglCreateContext(hdc);
                    if (!wglMakeCurrent(hdc, m_hrc)) {
                        wglDeleteContext(m_hrc);
                        m_hrc = NULL;
                    }
                    else {
                        int glewState = glewInit();
                        if (GLEW_OK != glewState) {
                            wglDeleteContext(m_hrc);
                            m_hrc = NULL;
                        }
                    }
                }
            }
            ReleaseDC(hWnd, hdc);
        }
	}

	WinOpenGLCanvas::~WinOpenGLCanvas() {
        if (m_hrc) {
            wglMakeCurrent(NULL, NULL);
            wglDeleteContext(m_hrc);
        }
	}

	bool WinOpenGLCanvas::IsValid() const { 
		return m_hrc; 
	}

	void WinOpenGLCanvas::Draw(Layout* layout) {
        if (m_hrc) {
            layout->Draw();
            HDC hdc = GetDC(m_hWnd);
            SwapBuffers(hdc);
            ReleaseDC(m_hWnd, hdc);
        }
	}

}

#endif
