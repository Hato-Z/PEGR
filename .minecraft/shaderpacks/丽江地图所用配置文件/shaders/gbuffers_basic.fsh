#version 120

varying vec4 color;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;


void main() {

	float skymask;
	float alphamask;
	
	if (gl_FragCoord.z < 0.992f) {
		skymask = 1.0f;
		alphamask = 0.0f;
	} else {
		skymask = 0.0f;
		alphamask = 1.0f;
	}
	

	gl_FragData[0] = color;
	gl_FragData[1] = vec4(vec3(gl_FragCoord.z), alphamask);
	


	
	
	gl_FragData[4] = vec4(0.0, 0.0, skymask, alphamask);
	
	//Illuminate black elements
	gl_FragData[5] = vec4(0.0, 0.0, 1.0, alphamask);
	gl_FragData[6] = vec4(0.0f, 0.0f, 0.0f, alphamask);

	if (fogMode == GL_EXP) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, (gl_Fog.color.rgb * 1.0), 1.0 - clamp(exp(-gl_Fog.density * gl_FogFragCoord), 0.0, 1.0));
	} else if (fogMode == GL_LINEAR) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, (gl_Fog.color.rgb * 1.0), clamp((gl_FogFragCoord - gl_Fog.start) * gl_Fog.scale, 0.0, 1.0));
	}
}