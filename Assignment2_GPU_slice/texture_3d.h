class Texture {

	public:
		Texture(GLenum TextureTarget){
			m_textureTarget = TextureTarget;
		}

	void Load3d(ScalarField *field ){
		
	
		glGenTextures(1, &m_textureObj);
		glBindTexture(m_textureTarget, m_textureObj);
		
		glTexImage3D(m_textureTarget, 0, GL_R32F, field->dim[0], field->dim[1],field->dim[2], 0, GL_RED, GL_FLOAT, field->data);
		
		glTexParameterf(m_textureTarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		
		glTexParameterf(m_textureTarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		
		glBindTexture(m_textureTarget, 0);
		
	}

	void Load1d(float *colorTexture){
	
		glGenTextures(1, &m_textureObj);
		glBindTexture(m_textureTarget, m_textureObj);
		
		glTexImage1D(m_textureTarget, 0, GL_RGB, 33, 0, GL_RGB, GL_FLOAT,colorTexture);
		
		glTexParameterf(m_textureTarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		
		glTexParameterf(m_textureTarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		
		glBindTexture(m_textureTarget, 0);
		
	}
	
	void Bind(GLenum TextureUnit)
	{
		glActiveTexture(TextureUnit);
		glBindTexture(m_textureTarget, m_textureObj);
	}

	private:
		GLenum m_textureTarget;
		GLuint m_textureObj;
};
