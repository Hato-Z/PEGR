#version 120

varying vec4 color;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;


void main() {

	gl_FragData[0] = color;
	gl_FragData[1] = vec4(vec3(gl_FragCoord.z), 1.0);
	
	float colorav = ((color.r + color.g + color.b) / 3.0);
	float alphamask = color.a;
	
	float colormask = 0.0;
	
		if (alphamask != 1.0 && colorav == 0.0) {
			colormask = 0.0;
		}  else {
			colormask = 1.0;
		}
		
		if (colorav == 1.0) {
			
		}
	
	gl_FragData[4] = vec4(0.0, 0.0, 1.0f - colormask, 1.0);

	if (fogMode == GL_EXP) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, (gl_Fog.color.rgb * 1.0), 1.0 - clamp(exp(-gl_Fog.density * gl_FogFragCoord), 0.0, 1.0));
	} else if (fogMode == GL_LINEAR) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, (gl_Fog.color.rgb * 1.0), clamp((gl_FogFragCoord - gl_Fog.start) * gl_Fog.scale, 0.0, 1.0));
	}
}