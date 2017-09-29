// Source.cpp : Defines the entry point for the console application.
//

#include <afxwin.h>  // necessary for MFC to work properly
#include "Header.h"
#include "../../src/blepo.h"
#include <stack>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace blepo;

// Calculate kernal widths
double calculateKernalHalfWidth(double sigma)
{
	return(floor(2.5*sigma) - 0.5);
}

double calculateKernalWidth(double sigma)
{
	return(floor(2 * (calculateKernalHalfWidth(sigma))) + 1);
}

// Calculate gaussians
double* gaussianKernal(double sigma)
{
	double width = calculateKernalWidth(sigma);
	double halfWidth = calculateKernalHalfWidth(sigma);
	double *Gaussian, sum = 0.0;

	Gaussian = new double[(int)width];
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = exp(-(i - halfWidth)*(i - halfWidth) / (2 * sigma*sigma));
		sum += Gaussian[i];
	}
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = (Gaussian[i] / sum);
	}

	return Gaussian;
	delete[] Gaussian;
}

double* derivativeKernal(double sigma)
{
	double width = calculateKernalWidth(sigma);
	double halfWidth = calculateKernalHalfWidth(sigma);
	double *Gaussian, sum = 0.0;

	Gaussian = new double[(int)width];
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = (i - halfWidth)*exp(-(i - halfWidth)*(i - halfWidth) / (2 * sigma*sigma));
		sum += (i*Gaussian[i]);
	}
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = (Gaussian[i] / sum);
	}

	return Gaussian;
	delete[] Gaussian;
}

// Convolution
ImgFloat convolve(ImgGray img_O, double* temp_1, double* temp_2, double sigma)
{
	ImgFloat img_1;
	img_1.Reset(img_O.Width(), img_O.Height());

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			img_1(x, y) = img_O(x, y);
		}
	}

	float sum;
	ImgFloat img_2, img_3;
	double width = calculateKernalWidth(sigma);
	double halfWidth = calculateKernalHalfWidth(sigma);

	img_2.Reset(img_O.Width(), img_O.Height());
	img_3.Reset(img_O.Width(), img_O.Height());

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			img_2(x, y) = 0;
			img_3(x, y) = 0;
		}
	}

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = halfWidth; x < img_O.Width() - halfWidth; ++x)
		{
			sum = 0;
			for (int z = 0; z < width; z++)
			{
				sum += ((temp_1[z] * img_1(x + halfWidth - z, y)));
			}
			img_2(x, y) = sum;
		}
	}

	for (int y = halfWidth; y < img_O.Height() - halfWidth; ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			sum = 0;
			for (int z = 0; z < width; z++)
			{
				sum += ((temp_2[z] * img_2(x, y + halfWidth - z)));
			}
			img_3(x, y) = sum;
		}
	}
	return img_3;
}

ImgFloat nonMaxSuppression(ImgFloat imgFile_1, ImgFloat imgFile_2)
{
	ImgFloat imgFile_output;
	double temp;
	imgFile_output.Reset(imgFile_1.Width(), imgFile_1.Height());
	for (int i = 0; i<imgFile_1.Width(); i++)
		for (int j = 0; j<imgFile_1.Height(); j++)
		{
			imgFile_output(i, j) = imgFile_1(i, j);
		}
	double pi = 3.14;
	for (int i = 1; i<imgFile_1.Width() - 1; i++)
	{
		for (int j = 1; j<imgFile_1.Height() - 1; j++)
		{
			temp = imgFile_2(i, j);
			if ((((-pi) / 8 <= temp&&temp<pi / 8) || ((-pi) / 8 <= temp + pi&&temp + pi<pi / 8)) && (imgFile_1(i, j)<imgFile_1(i - 1, j) && imgFile_1(i, j)<imgFile_1(i + 1, j)))
			{
				imgFile_output(i, j) = 0;
			}
			if ((((pi) / 8 <= temp&&temp<(3 * pi) / 8) || ((pi) / 8 <= temp + pi&&temp + pi<(3 * pi) / 8)) && (imgFile_1(i, j)<imgFile_1(i - 1, j - 1) && imgFile_1(i, j)<imgFile_1(i + 1, j + 1)))
			{
				imgFile_output(i, j) = 0;
			}
			if ((((3 * pi) / 8 <= temp&&temp<(5 * pi) / 8) || ((3 * pi) / 8 <= temp + pi&&temp + pi<(5 * pi) / 8)) && (imgFile_1(i, j)<imgFile_1(i, j - 1) && imgFile_1(i, j)<imgFile_1(i, j + 1)))
			{
				imgFile_output(i, j) = 0;
			}
			if ((((5 * pi) / 8 <= temp&&temp<(7 * pi) / 8) || ((5 * pi) / 8 <= temp + pi&&temp + pi<(7 * pi) / 8)) && (imgFile_1(i, j)<imgFile_1(i - 1, j + 1) && imgFile_1(i, j)<imgFile_1(i + 1, j - 1)))
			{
				imgFile_output(i, j) = 0;
			}
		}
	}
	return imgFile_output;
}


ImgBinary performThreshold(ImgFloat img_1, int temp)
{
	ImgBinary img_2;
	img_2.Reset(img_1.Width(), img_1.Height());

	for (int y = 0; y < img_1.Height(); ++y)
	{
		for (int x = 0; x < img_1.Width(); ++x)
		{
			img_2(x, y) = img_1(x, y);
			if (img_1(x, y) > temp)
			{
				img_2(x, y) = 1;
			}
			else
			{
				img_2(x, y) = 0;
			}
		}
	}
	return img_2;
}

void chamferImage(ImgBinary img_1, ImgInt *img_2)
{
	img_2->Reset(img_1.Width(), img_1.Height());
	int temp_1 = (img_1.Width()*img_1.Height()) + 1;

	for (int y = 0; y < img_1.Height(); ++y)
	{
		for (int x = 0; x < img_1.Width(); ++x)
		{
			if (img_1(x, y))
			{
				(*img_2)(x, y) = 0;
			}
			else
			{
				int temp_2 = temp_1;
				if (y > 0)
					temp_2 = blepo_ex::Min(temp_2, (*img_2)(x, y - 1) + 1);
				if (x > 0)
					temp_2 = blepo_ex::Min(temp_2, (*img_2)(x - 1, y) + 1);
				(*img_2)(x, y) = temp_2;
			}
		}
	}

	for (int y = img_1.Height() - 1; y >= 0; y--)
	{
		for (int x = img_1.Width() - 1; x >= 0; x--)
		{
			if (img_1(x, y))
				(*img_2)(x, y) = 0;
			else
			{
				int temp_3 = (*img_2)(x, y);
				if (y < img_1.Height() - 1)
					temp_3 = blepo_ex::Min(temp_3, (*img_2)(x, y + 1) + 1);
				if (x < img_1.Width() - 1)
					temp_3 = blepo_ex::Min(temp_3, (*img_2)(x + 1, y) + 1);
				(*img_2)(x, y) = temp_3;
			}
		}
	}
}

void floodfillImage(ImgBinary img_1, int a, int b, int z, ImgBinary* img_2)
{
	int temp = img_1(a, b);
	std::stack<Point> temp_1;
	Point temp_2, seedP;
	temp_2.SetPoint(a, b);

	if (temp == z)
		return;

	temp_1.push(Point(a, b));
	(*img_2)(a, b) = z;

	while (temp_1.empty())
	{
		seedP = temp_1.top();
		int x = seedP.x;
		int y = seedP.y;
		temp_1.pop();

		if (img_1(x + 1, y) == temp && (*img_2)(x + 1, y) != z)
		{
			temp_1.push(Point(x + 1, y));
			(*img_2)(x + 1, y) = z;
		}
		if (img_1(x - 1, y) == temp && (*img_2)(x - 1, y) != z)
		{
			temp_1.push(Point(x - 1, y));
			(*img_2)(x - 1, y) = z;
		}
		if (img_1(x + 1, y - 1) == temp && (*img_2)(x + 1, y - 1) != z)
		{
			temp_1.push(Point(x + 1, y - 1));
			(*img_2)(x + 1, y - 1) = z;
		}
		if (img_1(x - 1, y + 1) == temp && (*img_2)(x - 1, y + 1) != z)
		{
			temp_1.push(Point(x - 1, y + 1));
			(*img_2)(x - 1, y + 1) = z;
		}
		if (img_1(x, y + 1) == temp && (*img_2)(x, y + 1) != z)
		{
			temp_1.push(Point(x, y + 1));
			(*img_2)(x, y + 1) = z;
		}
		if (img_1(x, y - 1) == temp && (*img_2)(x, y - 1) != z)
		{
			temp_1.push(Point(x, y - 1));
			(*img_2)(x, y - 1) = z;
		}
		if (img_1(x + 1, y + 1) == temp && (*img_2)(x + 1, y + 1) != z)
		{
			temp_1.push(Point(x + 1, y + 1));
			(*img_2)(x + 1, y + 1) = z;
		}
		if (img_1(x - 1, y - 1) == temp && (*img_2)(x - 1, y - 1) != z)
		{
			temp_1.push(Point(x - 1, y - 1));
			(*img_2)(x - 1, y - 1) = z;
		}
	}
	return;
}



ImgBinary performDoubleThreshold(ImgFloat imgFile_input, double low, double high)
{
	ImgBinary imgFile_output, imgFile_output1, imgFile_tLow;
	imgFile_output1.Reset(imgFile_input.Width(), imgFile_input.Height());
	imgFile_output.Reset(imgFile_input.Width(), imgFile_input.Height());
	imgFile_tLow.Reset(imgFile_input.Width(), imgFile_input.Height());

	imgFile_tLow = performThreshold(imgFile_input, low);

	for (int x = 0; x<imgFile_input.Width(); x++)
	{
		for (int y = 0; y<imgFile_input.Height(); y++)
		{
			if (imgFile_output1(x, y) == 1)
			{
				imgFile_output1(x, y) = 0; imgFile_output(x, y) = 0;
			}
		}
	}

	for (int x = 0; x<imgFile_input.Width(); x++)
	{
		for (int y = 0; y<imgFile_input.Height(); y++)
		{
			if (imgFile_input(x, y)>low)
			{
				imgFile_tLow(x, y) = 1;
				imgFile_output1(x, y) = 0;
			}
		}
	}

	for (int x = 0; x<imgFile_input.Width(); x++)
	{
		for (int y = 0; y<imgFile_input.Height(); y++)
		{
			if (imgFile_input(x, y)>high)
			{
				floodfillImage(imgFile_output, x, y, 1, &imgFile_output1);
			}
		}
	}

	ImgInt chamfer;
	chamfer.Reset(imgFile_output1.Width(), imgFile_output1.Height());

	return imgFile_output1;
}

ImgBinary performEdgeLink(ImgFloat imgFile_img1)
{
	int i, j, count = 0, k = 0;
	for (i = 0; i<imgFile_img1.Width(); i++)
		for (j = 0; j<imgFile_img1.Height(); j++)
			count++;
	double *array;
	array = new double[(int)count];
	for (i = 0; i<imgFile_img1.Width(); i++)
	{
		for (j = 0; j<imgFile_img1.Height(); j++)
		{
			if (k <= count)
			{
				array[k] = imgFile_img1(i, j);
				k++;
			}
		}
	}

	std::sort(array, array + count, std::greater<int>());
	int *temp, temp_1, size;
	size = floor(array[0]);
	temp = new int[(int)size];
	for (i = 0; i <= size; i++)
		temp[i] = 0;
	for (i = 0; i<count; i++)
	{
		if (array[i] != NULL)
		{
			temp_1 = floor(array[i]);
			temp[temp_1]++;
		}
	}
	int *cdf, sum = 0;
	cdf = new int[(int)size];
	for (i = size; i >= 0; i--)
	{
		sum = sum + temp[i];
		cdf[i] = sum;
	}
	int temp_2 = 10 * cdf[0] / 100;
	double highT = 0;
	for (i = 0; i <= size; i++)
	{
		if (temp_2 + 500>cdf[i] && temp_2 - 500<cdf[i])
		{
			highT = i;
		}
	}
	double LowT = highT / 5;
	ImgBinary imgFile_out = performDoubleThreshold(imgFile_img1, LowT, highT);

	delete array;
	delete temp;
	delete cdf;
	return imgFile_out;

}

void ImgGradient(ImgGray img_O, double sigma, ImgFloat *img_1, ImgFloat *img_2, ImgFloat *img_3, ImgFloat *img_4)
{
	(*img_1).Reset(img_O.Width(), img_O.Height());
	(*img_2).Reset(img_O.Width(), img_O.Height());

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			(*img_1)(x, y) = 0;
			(*img_2)(x, y) = 0;
		}
	}

	double *Gaussian = gaussianKernal(sigma);
	double *DGaussian = derivativeKernal(sigma);
	double *flip;
	double width = calculateKernalWidth(sigma);

	flip = new double[(int)width];
	int temp = width - 1;
	for (int z = 0; z < width; z++)
	{
		flip[temp] = DGaussian[z];
		temp--;
	}

	(*img_3) = convolve(img_O, flip, Gaussian, sigma);
	(*img_4) = convolve(img_O, Gaussian, flip, sigma);


	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			(*img_1)(x, y) = sqrt(((*img_3)(x, y)*(*img_3)(x, y)) + (((*img_4)(x, y)*(*img_4)(x, y))));
			(*img_2)(x, y) = atan2((*img_4)(x, y), (*img_3)(x, y));
		}
	}
	return;
}


int main(int argc, const char* argv[], const char* envp[])
{
	// Initialize MFC and return if failure
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", (int)hModule);
		return 1;
	}

	if (argc<3 || argc>4)
	{
		cout << "Error! argument count is invalid" << endl;
		exit(0);
	}

	try
	{
		double sigma;
		sigma = atof(argv[1]);

		ImgGray imgFile_O1;
		ImgBgr imgFile_O2;

		CString path("../../images/");
		CString imgPath_1 = path + CString(argv[2]);

		Load(imgPath_1, &imgFile_O1);
		Load(imgPath_1, &imgFile_O2);

		Figure fig;
		fig.SetTitle("Input Image");
		fig.Draw(imgFile_O2);

		double width = calculateKernalWidth(sigma);
		double *Gaussian = gaussianKernal(sigma);
		double *DGaussian = derivativeKernal(sigma);

		cout << "The Guassian Kernel is ";
		for (int i = 0; i<width; i++)
		{
			cout << "\t" << Gaussian[i];
		}

		cout << endl << "The Derivative of the Guassian Kernel is";
		for (int i = 0; i<width; i++)
		{
			cout << "\t" << DGaussian[i];
		}

		ImgFloat imgFile_mag, imgFile_phase, imgFile_x, imgFile_y;
		ImgGradient(imgFile_O1, sigma, &imgFile_mag, &imgFile_phase, &imgFile_x, &imgFile_y);

		Figure grad_x, grad_y;
		grad_x.SetTitle("X-Gradient");
		grad_y.SetTitle("Y-Gradient");
		grad_x.Draw(imgFile_x);
		grad_y.Draw(imgFile_y);

		Figure mag, phase;
		mag.SetTitle("Magnitude Gradient");
		mag.Draw(imgFile_mag);
		phase.SetTitle("Phase Gradient");
		phase.Draw(imgFile_phase);

		ImgFloat imgFile_nonMaxSup = nonMaxSuppression(imgFile_mag, imgFile_phase);
		Figure nonMaxSup;
		nonMaxSup.SetTitle("Non Maximum Suppression");
		nonMaxSup.Draw(imgFile_nonMaxSup);

		ImgBinary imgFile_edge = performEdgeLink(imgFile_nonMaxSup);
		Figure edge;
		edge.SetTitle("Edge Linked Image");
		edge.Draw(imgFile_edge);

		ImgInt imgFile_int;
		chamferImage(imgFile_edge, &imgFile_int);
		Figure chamf;
		chamf.SetTitle("Chamfered Image");
		chamf.Draw(imgFile_int);

		EventLoop();
	}

	catch (Exception& e)
	{
		e.Display();
	}

	return 0;
}
