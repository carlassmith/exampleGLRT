import ij.*;
import ij.IJ;
import ij.gui.*;

import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.frame.Fitter;
import ij.plugin.Duplicator;

import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.FHT;

import ij.measure.CurveFitter;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;

import java.lang.Math;
import java.util.Arrays;
import java.awt.*;

/*
 
(C) Copyright 2016
All rights reserved             Faculty of Applied Sciences
                                Delft University of Technology
                                Delft, The Netherlands
 
 Bernd Rieger & Diederik Feilzer
 email: b.rieger@tudelft.nl
 
 Literature: 
    R. Heintzmann, P. Relich, R. Nieuwenhuizen, K. Lidke, B. Rieger,
    Calibrating photon counts from a single image, submitted
 
 */

public class Image_pcfo implements PlugIn {

	// image property members
	private int input_width;
	private int input_height;

	// plugin parameters
	public float kcut;
	public float RNStd;
	public boolean AvoidCross;
	
	//fixed params
	float borderfraction = 0.05f; 
	int CrossWidth = 3; 
 
    public void run(String arg) {
    
    	if (arg.equals("about")) {
			IJ.showMessage("About", "Calculate gain and offset from a single image.");
			return;
		}
    
        //get the active image
        ImagePlus impin = IJ.getImage();
        ImagePlus imp = new Duplicator().run(impin);
                
        //convert to grey 32-bit float if not all ready
        ImageConverter ic = new ImageConverter(imp);
		ic.convertToGray32();
        
        if (null == imp)
        	return;
        
        //ask for settings
        GenericDialog gd = new GenericDialog("Parameters");
        gd.addNumericField("kcut:", 0.9, 2);
        gd.addNumericField("RNStd:", 0.0, 2);
        gd.addNumericField("AvoidCross:", 1, 0);
        gd.addNumericField("TestMode:", 0, 0);
        gd.showDialog();
        
        if (gd.wasCanceled())
        	return;
 
        kcut = (float)gd.getNextNumber();
        RNStd = (float)gd.getNextNumber();
        AvoidCross = ((int)gd.getNextNumber()==1);
		
		ImageProcessor ip = imp.getStack().getProcessor(1);
		
		if ((int)gd.getNextNumber()==1) {
			ip = addNoiseAndOffset(ip);
		}
		
		input_width = ip.getWidth();
		input_height = ip.getHeight();
		
		//start main function
		ImageProcessor p = pcfo_TileOnly(ip);
		
    }
    
	private ImageProcessor addNoiseAndOffset(ImageProcessor ip) {
   
   		GenericDialog gd = new GenericDialog("Parameters");
        gd.addNumericField("gain:", 0.5, 2);
        gd.addNumericField("offset:", 0.0, 1);
        gd.addNumericField("PSF width:", 2.0, 1);
        gd.showDialog();
        
        if (gd.wasCanceled())
        	return ip;
        	
        float gain = (float)gd.getNextNumber();
        float offset = (float)gd.getNextNumber();
        double psfw = (double)gd.getNextNumber();
   
		int width = ip.getWidth(); // width of the image
		int height = ip.getHeight(); // height of the image
		
		// simple simulation of PSF.
		GaussianBlur gs = new GaussianBlur();
		gs.blurGaussian(ip,(double)(psfw), (double)(psfw), 0.02);

		float[] pixels = (float[]) ip.getPixels();

		for (int i = 0; i < width*height; i++) {
			pixels[i] = poissonrnd(pixels[i]*gain)/gain + offset;
		}
          
		return ip;
	}
	
	/*
	https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
	For mu > 5 it takes very long.
	*/

	private int poissonrnd(float mu) {

	double L = Math.exp(-mu);
	int k = -1;
	double p = 1;
		while (p >= L) {
			k++;
			p *= Math.random();
		}
		return k;
	}
    
    private ImageProcessor pcfo_TileOnly(ImageProcessor ip) {
    	
    	int TilesX = 3;//number of tiles, hard coded - in matlab code this is a parameter
		int TilesY = 3;
		int n = 0;
		
		int width = ip.getWidth();
		int height = ip.getHeight();
		
		int[] 	NumPixelsInBin 	= new int[TilesX * TilesY];
		float[] TotalVar 		= new float[TilesX * TilesY];
		float[] TotalInt 		= new float[TilesX * TilesY];

		//loop trough tiles
		
		for (int ty=0; ty < TilesY; ty ++) {
    		for (int tx=0; tx < TilesX; tx ++) {
        		int c_x 	= (int)Math.floor((tx*width)/TilesX);
        		int c_w 	= (int)Math.floor(width/TilesX);
        		int c_y 	= (int)Math.floor((ty*height)/TilesY);
        		int c_h 	= (int)Math.floor(height/TilesY);
        
        		if(c_h*c_w > 0) {
        			NumPixelsInBin[n] = c_h*c_w;
        			ImageProcessor tile = crop(ip, c_x, c_y, c_w, c_h);

					float[] p = pcf(tile);
        			
        			TotalVar[n] = p[0];
        			TotalInt[n] = p[1];
        			
        			n++;
        		}     
        		
        	}
        }
        
        //manage cases of missing values
        TotalVar = Arrays.copyOfRange(TotalVar, 1, n);
        TotalInt = Arrays.copyOfRange(TotalInt, 1, n);
        NumPixelsInBin = Arrays.copyOfRange(NumPixelsInBin, 1, n);
        
        //fit a straight line
        float[] respara = RWLSPoisson(TotalInt,TotalVar,NumPixelsInBin);
       
		float Voffset = respara[0];
		float slope = respara[1];
		
		//take pcfRH over whole image to avoid tile effects.
		float[] po = pcf(ip);
		
		float AllTotalVar = po[0];
		float AllTotalInt = po[1];
		
		float gain = ((float)(AllTotalInt))/(((float)(AllTotalVar)-Voffset));
		
		float offset = (RNStd*RNStd/gain)-Voffset/slope;
		
		IJ.showMessage("results","Offset: " + String.format("%.2f", (float)offset) + "\n Gain: " + String.format("%.3f", 1/(float)gain));
				
		return ip;
    }
    
    /*
	Effectively places a mirror on the right and bottom edge of the image.
	*/
    
    private ImageProcessor symmetrize(ImageProcessor ip) {
        
        int width = ip.getWidth();
		int height = ip.getHeight();
        
        float[] pixels_in = (float[]) ip.getPixels();
        
        float[] pixels_out = new float[2*width * 2*height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
				float pixel_in = pixels_in[x + y * width];
                pixels_out[x+2*width*y] 							= pixel_in;
                pixels_out[2*width-x-1+2*width*y] 					= pixel_in;
                pixels_out[x+2*width*(2*height-y-1)]				= pixel_in;
                pixels_out[(2*width-x-1)+2*width*(2*height-y-1)] 	= pixel_in;
            }
        }
        return new FloatProcessor(2*width, 2*height, pixels_out, null);
    }
    
    /*
	Crops ImageProcessor to specified dimensions
	*/
    
    private ImageProcessor crop(ImageProcessor ip, int c_x, int c_y, int c_w, int c_h) {
        
        int width = ip.getWidth();
		int height = ip.getHeight();
        
        float[] pixels_in = (float[]) ip.getPixels();
        
        float[] pixels_out = new float[c_w * c_h];
        for (int y = c_y; y < c_y+c_h; y++) {
            for (int x = c_x; x < c_x+c_w; x++) {

                pixels_out[x-c_x+c_w*(y-c_y)] = pixels_in[x + y * width];
   
            }
        }
        
        return new FloatProcessor(c_w, c_h, pixels_out, null);
    }
    
    /*
	...
	*/
    
    private float[] pcf(ImageProcessor ip) {

		ImageProcessor ip_sym = symmetrize(ip);
		
		ImageProcessor ip_filt = ft_mask(ip_sym);
		
		ImageProcessor amask = mask2(ip);
		
		ImageProcessor masker = mask(ip_filt);
				
		float fra = mean(masker);

		float TotalVar = mean(ip_filt,masker)/fra;
		float ImgSum = mean(ip,amask);
		
    	return new float[] {(float)TotalVar, (float)ImgSum};
    }
    
    /*
	...
	*/
    
    private ImageProcessor ft_mask(ImageProcessor ip) {
    
    	int width = ip.getWidth();
		int height = ip.getHeight();
    
        int maxN = Math.max(width, height);
        int size = 2;
        while(size<maxN){
        	size *= 2;
        }

        int h = Math.round((size-height)/2);
        int w = Math.round((size-width)/2);
        
        // create padded version of ip in ip2
        ImageProcessor ip2 = zeros(size, size);  
        ip2.insert(ip, w, h);
        
        offset(ip2,-mean(ip2));
        
        // calculate the Fourier transformation
    	FHT fht = new FHT(ip2);
        fht.transform();	
        
        ImageStack amp = fht.getComplexTransform();
                
        float[] pixels_real = (float[])amp.getProcessor(1).getPixels();
        float[] pixels_imag = (float[])amp.getProcessor(2).getPixels();
        
      	float[] pixels_sqr = new float[size*size];
      	
        for (int i=0; i<size*size; i++) {
    	    pixels_sqr[i] = (float)((pixels_real[i]*pixels_real[i]+pixels_imag[i]*pixels_imag[i])/(width*height*Math.sqrt(2*Math.PI)));
        }
        
        FloatProcessor ampsqr = new FloatProcessor(size, size, pixels_sqr, null);
        
        return ampsqr;
    }
    
	/*
	...
	*/
    
    private ImageProcessor mask(ImageProcessor ip) {
    
    	int width = ip.getWidth();
		int height = ip.getHeight();
    
    	float radius = (float)(kcut*height/2.0);
				
		float[] mask_pixels = new float[width * height];
		
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                mask_pixels[x+width*y] = 1f;
                
                float rx = (float)x-width/2;
				float ry = (float)y-height/2;
                
                if(Math.sqrt(rx*rx+ry*ry) < radius) {
    				mask_pixels[x+width*y] = 0f;
    			}
    			                
                if (AvoidCross) {
    				if((Math.abs(rx) <= CrossWidth) || (Math.abs(ry) <= CrossWidth)) {
    					mask_pixels[x+width*y] = 0f;
    				}
				}
				
            }
        }
        
        FloatProcessor fp = new FloatProcessor(width, height, mask_pixels, null);
        
        GaussianBlur gs = new GaussianBlur();
        
        gs.blurGaussian(fp,(double)(CrossWidth/4), (double)(CrossWidth/4), 0.01);
        
        return fp;
    }
    
    /*
	...
	*/
    
    private ImageProcessor mask2(ImageProcessor ip) {
    
    	int width = ip.getWidth();
		int height = ip.getHeight();
    
    	int border = (int)(borderfraction/2*input_height);
				
		float[] mask_pixels = new float[width * height];
		
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                mask_pixels[x+width*y] = 1f;
                
                if((x <= border) || (y <= border) || (x >= width-border) || (y >= height-border)) {
    				mask_pixels[x+width*y] = 0f;
    			}
    
            }
        }
        
        FloatProcessor fp = new FloatProcessor(width, height, mask_pixels, null);

        return fp;
    }
    
    /*
	Returns an zero filled ImageProcessor.
	*/
    
    private ImageProcessor zeros(int width, int height) {
    
		float[] pixels = new float[width * height];
		
        for (int i = 0; i < height * width; i++) {
            pixels[i] = 0f;
        }
        
        FloatProcessor fp = new FloatProcessor(width, height, pixels, null);

        return fp;
    }
    
    private ImageProcessor abssqr(ImageProcessor ip) {
    	
    	int width = ip.getWidth();
		int height = ip.getHeight();
        
        float[] pixels_in = (float[]) ip.getPixels();
                        
        for (int i = 0; i < height * width; i++) {
			pixels_in[i] = pixels_in[i]*pixels_in[i];
        }
    
    	return new FloatProcessor(width, height, pixels_in, null);
    }
    
    private float sum(ImageProcessor ip) {
    	
    	int width = ip.getWidth();
		int height = ip.getHeight();
        
        float[] pixels_in = (float[]) ip.getPixels();
        
        float sum = 0f;
        
        for (int i = 0; i < height * width; i++) {
			sum += pixels_in[i];
        }
    
    	return (float) sum;
    }
    
	private float mean(ImageProcessor ip) {
    	
    	int width = ip.getWidth();
		int height = ip.getHeight();
    
    	return (float) sum(ip)/(width*height);
    }
    
    private float mean(ImageProcessor ip, ImageProcessor mask) {
    	
    	int width = ip.getWidth();
		int height = ip.getHeight();
		
		int m_width = mask.getWidth();
		int m_height = mask.getHeight();
		
		if ((m_width != width) || (m_height != height))
			return 0f;
        
        float[] pixels_in = (float[]) ip.getPixels();
        float[] mask_in = (float[]) mask.getPixels();
        
        float sum = 0f;
        float m_sum = 0f;
        
        for (int i = 0; i < height * width; i++) {
        	float m = mask_in[i];
			sum += pixels_in[i]*m;
			m_sum += m;
        }
    
    	return (float) sum/m_sum;
    }
    
    private ImageProcessor offset(ImageProcessor ip, float offset) {
    	
    	int width = ip.getWidth();
		int height = ip.getHeight();
        
        float[] pixels_in = (float[]) ip.getPixels();
                
        for (int i = 0; i < height * width; i++) {
			pixels_in[i] = pixels_in[i]+offset;
        }
    
    	return new FloatProcessor(width, height, pixels_in, null);
    }
    
    /*
    Poisson Reweighted Least Square fit, linear regression with error according to the Poisson distribution.
    Fits a straight line with offset.
    */
    
    private float[] RWLSPoisson(float[] x, float[] y, int[] N) {
    	
    	int l = x.length;
    	float o = 0f;
    	float s = 0f;
    	float myTresh = 1f;
    	int NumIter = 5;
    	float[] v = new float[l];
    	float[] vv = new float[l];
    	v = y;
    	
    	for (int n = 0; n < NumIter; n++) {
        	for (int j = 0; j < l; j++) {
        		vv[j] = v[j]*v[j]/N[j];
        	}
        	boolean er = false;
        	for (int j = 0; j < l; j++) {
        		if (v[j] < myTresh){
        			er = true;
        			v[j] = myTresh;
        		}
        	}
        	if (er) {
        		if (n==1) {
        			IJ.showMessage("Warning:","The data has a variance below 2 ADUs at low signal level. This leads to unwanted biases.\nIncreasing the variance estimation for the fit.");
        		}
        	}
        	float[] p = GLS(x,y,v);
        	o = p[0];
        	s = p[1];
        	
        	for (int j = 0; j < l; j++) {
        		v[j] = o + s*x[j];
        	}
        	
        }
        
        return new float[] {(float)o, (float)s};
    	
    }
    
    /*
    Generalized Least Square fit, linear regression with error.
    Fits a straight line with offset through data with known variances.
    */
    
    private float[] GLS(float[] x, float[] y, float[] v) {	
    	
    	int l = x.length;
    	float a = 0f;
    	float b = 0f;
    	float d = 0f;
    	float e = 0f;
    	float f = 0f;
        
        for (int i = 0; i < l; i++) {
            a+=1/v[i];
            b+=x[i]/v[i];
            d+=x[i]*x[i]/v[i];
           	e+=y[i]/v[i];
            f+=x[i]*y[i]/v[i];
        }
        
        float k = a*d-b*b;
        float o = (d*e-b*f) / k;
        float s = (a*f-b*e) / k;
    
    	return new float[] {(float)o, (float)s};
    }
 
}
