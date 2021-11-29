#include <omp.h>
#include <endian.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#define calcIndex(width, x, y) ((y) * (width) + (x))

int countNeighbors(double *currentfield, int x, int y, int w, int h);

long TimeSteps = 3;

void writeVTK2(long timestep, double *data, char prefix[1024], int localWidth, int localHeight, int threadNumber, int originX, int originY, int totalWidth)
{
    char filename[2048];
    int x, y;
    int w = localWidth - originX, h = localHeight - originY;
    int offsetX = originX;
    int offsetY = originY;
    float deltax = 1.0;
    long nxy = w * h * sizeof(float);
    //printf("NXY:%ld\n", nxy);
    snprintf(filename, sizeof(filename), "%s%s_%d-%05ld%s", "vtiFiles/", prefix, threadNumber, timestep, ".vti");
    FILE *fp = fopen(filename, "w");

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(fp, "<ImageData WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX,
            offsetX + w, offsetY, offsetY + h, 0, 0, deltax, deltax, 0.0);
    fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
    fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
    fprintf(fp, "</CellData>\n");
    fprintf(fp, "</ImageData>\n");
    fprintf(fp, "<AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");
    fwrite((unsigned char *)&nxy, sizeof(long), 1, fp);

    // printf("\nThread: %d writing...\noriginX/originY:%d/%d localWidth/localHeight: %d/%d\n", threadNumber, originX, originY, w, h);

    for (y = originY; y < localHeight; y++)
    {
        for (x = originX; x < localWidth; x++)
        {
            float value = data[calcIndex(totalWidth, x, y)];
            // if (threadNumber == 0)
            //     printf("Entry: i: %d, %f | thread:%d\n", calcIndex(totalWidth, x, y), value, threadNumber);
            fwrite((unsigned char *)&value, sizeof(float), 1, fp);
        }
    }
    //printf("\nThread: %d finished\n", threadNumber);

    fprintf(fp, "\n</AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

void show(double *currentfield, int w, int h)
{
    printf("\033[H");
    int x, y;
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
            printf(currentfield[calcIndex(w, x, y)] ? "\033[07m  \033[m" : "  ");
        //printf("\033[E");
        printf("\n");
    }
    fflush(stdout);
}

int countNeighbors(double *currentfield, int x, int y, int widthTotal, int heightTotal)
{
    int cnt = 0;
    int y1 = y - 1;
    int x1 = x - 1;

    while (y1 <= y + 1)
    {
        if (currentfield[calcIndex(widthTotal, (x1 + widthTotal) % widthTotal, (y1 + heightTotal) % heightTotal)])
            cnt++;

        if (x1 <= x)
            x1++;
        else
        {
            y1++;
            x1 = x - 1;
        }
    }
    return cnt;
}

void evolve(double *currentfield, double *newfield, int widthTotal, int heightTotal, int threadNumber, int nX, int nY, int pX, int pY, long t)
{
    int x, y, xStart;
    int localpX, localpY;

    localpX = threadNumber % pX;
    localpY = threadNumber / pX;

    x = nX * localpX;
    xStart = x;

    int widthLocal = nX * (localpX + 1);
    y = nY * localpY;
    int heightLocal = nY * (localpY + 1);

    int originX = x, originY = y;

    //printf("Thread: %d, localpX: %d, localpY: %d, X:%d,widthLocal: %d, Y:%d, heightLocal: %d\n",threadNumber,localpX,localpY,x,widthLocal,y,heightLocal);
    while (y < heightLocal)
    {
        int neighbors = countNeighbors(currentfield, x, y, widthTotal, heightTotal);
        int index = calcIndex(widthTotal, x, y);
        if (currentfield[index])
            neighbors--;
        //if (neighbors >= 2)
        //printf("widthTotal: %d | x/y: %d/%d | Index:%d Neighbors: %d Thread: %d, Am I? %f\n", widthTotal, x, y, index, neighbors, threadNumber, currentfield[index]);

        if (neighbors == 3 || (neighbors == 2 && currentfield[index]))
        {
            newfield[index] = 1;
            // printf("Test: %f, X/Y: %d/%d\n", newfield[calcIndex(widthTotal, x, y)],x,y);
        }
        else
        {
            newfield[index] = 0;
        }
        if (x < widthLocal - 1)
            x++;
        else
        {
            y++;
            x = xStart;
        }
    }
    writeVTK2(t, currentfield, "gol", widthLocal, heightLocal, threadNumber, originX, originY, widthTotal);
}

void filling(double *currentfield, int w, int h)
{
    int i;

    for (i = 0; i < h * w; i++)
    {
        currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
}

void game(int nX, int nY, int threadX, int threadY)
{
    int widthTotal = nX * threadX;
    int heightTotal = nY * threadY;

    double *currentfield = calloc(widthTotal * heightTotal, sizeof(double));
    double *newfield = calloc(widthTotal * heightTotal, sizeof(double));

    filling(currentfield, widthTotal, heightTotal);

    long t;
    for (t = 0; t < TimeSteps; t++)
    {
        // show(currentfield, widthTotal, heightTotal);

#pragma omp parallel num_threads(threadX *threadY)
        {
            int threadNumber = omp_get_thread_num();
            //printf("Max threads: %d | Used threads: %d | ThreadNr.: %d\n", omp_get_max_threads(), threadY * threadX, threadNumber);
            evolve(currentfield, newfield, widthTotal, heightTotal, threadNumber, nX, nY, threadX, threadY, t);
        }

        //printf("%ld timestep\n", t);

        //writeVTK2(t, currentfield, "gol", widthTotal, heightTotal);

        //usleep(100000);

        ///SWAP
        double *temp = currentfield;
        currentfield = newfield;
        newfield = temp;
    }

    free(currentfield);
    free(newfield);
}

int main(int argc, char *argv[])
{
    int segmentWidth = 0;
    int segmentHeight = 0;
    int amountXThreads = 0;
    int amountYThreads = 0;

    if (argc > 0)
        TimeSteps = atoi(argv[1]);
    if (argc > 1)
        segmentWidth = atoi(argv[2]);
    if (argc > 2)
        segmentHeight = atoi(argv[3]);
    if (argc > 3)
        amountXThreads = atoi(argv[4]);
    if (argc > 4)
        amountYThreads = atoi(argv[5]);

    // default:
    if (segmentWidth <= 0)
        segmentWidth = 200;
    if (segmentHeight <= 0)
        segmentHeight = 200;
    if (amountXThreads <= 0)
        amountXThreads = 2;
    if (amountYThreads <= 0)
        amountYThreads = 4;
    if (TimeSteps <= 0)
        TimeSteps = 100;

    game(segmentWidth, segmentHeight, amountXThreads, amountYThreads);
}
