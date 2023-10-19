#version 120

uniform sampler2D texture;

varying vec4 color;
varying vec4 texcoord;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;

void main() {

	float fullalpha = (texture2D(texture, texcoord.st).a * color.a);

	gl_FragData[0] = vec4(texture2D(texture, texcoord.st).rgb * color.rgb, fullalpha*1.0f);
	gl_FragData[1] = vec4(vec3(gl_FragCoord.z*1.0f), 1.0);
	
	float colormask = 0.0;
	float coloraverage = (color.r + color.g + color.b)/3.0;
	
	if (coloraverage == 1.0 && gl_FragCoord.z < 0.999) {
		colormask = 1.0;
	} else {
		colormask = 0.0;
	}
	
	gl_FragData[4] = vec4(0.0, 0.0, colormask, 1.0);
	gl_FragData[5] = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	//gl_FragData[1] = vec4(0.0);
		
	if (fogMode == GL_EXP) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, gl_Fog.color.rgb, 1.0 - clamp(exp(-gl_Fog.density * gl_FogFragCoord), 0.0, 1.0));
	} else if (fogMode == GL_LINEAR) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, gl_Fog.color.rgb, clamp((gl_FogFragCoord - gl_Fog.start) * gl_Fog.scale, 0.0, 1.0));
	}
}