import edu.princeton.cs.algs4.Picture;
import java.awt.Color;

public class SeamCarver {
    private Picture pic;
    private double [][] energyGrid;
    private int width, height;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        pic = new Picture(picture);
        width = pic.width();
        height = pic.height();
        energyGrid = new double[width][height];

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                energyGrid[i][j] = energyCal(i, j);
            }
        }
        
    }  

    // current picture
    public Picture picture() {
        Picture picture = new Picture(width, height);
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                picture.set(i, j, pic.get(i, j));
            }
        }
        return picture;
    }

    // width of current picture
    public int width() {
        return width;
    }

    // height of current picture
    public int height() {
        return height;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        return energyGrid[x][y];
    }

    private double energyCal(int x, int y) {
        if (border(x, y)) return 1000;
        return Math.sqrt(colorDiff(x - 1, y, x + 1, y) 
               + colorDiff(x, y - 1, x, y + 1));
    }

    private boolean border(int x, int y) {
        if ((x == 0) || (y == 0)) return true;
        if ((x == width - 1) || (y == height - 1)) return true;
        return false;
    }

    private double colorDiff(int px, int py, int nx, int ny) {
        Color p = pic.get(px, py);
        Color n = pic.get(nx, ny);
        int redDiff = p.getRed() - n.getRed();
        int greenDiff = p.getGreen() - n.getGreen();
        int blueDiff = p.getBlue() - n.getBlue();
        return redDiff*redDiff + greenDiff*greenDiff + blueDiff*blueDiff;
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        double[] distTo = new double[height];
        int[][] edgeTo = new int[width][height];
        for (int j = 0; j < width; j++) {
            double[] newDist = new double[height];
            for (int i = 0; i < height; i++) {
                if (j == 0) {
                    distTo[i] = energyGrid[j][i];
                } else {
                    double left, middle, right;
                    if (i == 0) left = 1000*width;
                    else left = distTo[i - 1];
                    middle = distTo[i];
                    if (i == height - 1) right = 1000*width;
                    else right = distTo[i + 1];
                    edgeTo[j][i] = i + min(left, middle, right);
                    newDist[i] = energyGrid[j][i] + distTo[edgeTo[j][i]];
                }
            }
            distTo = newDist;
        }
        double minDist = 1000*width;
        int[] minRoute = new int[width];
        for (int i = 0; i < height; i++) {
            if (distTo[i] < minDist) {
                minDist = distTo[i];
                minRoute[width - 1] = i;
            }
        }

        for (int j = width - 2; j >= 0; j--) {
            minRoute[j] = edgeTo[j + 1][minRoute[j + 1]];
        }

        return minRoute;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        double[] distTo = new double[width];
        int[][] edgeTo = new int[width][height];
        for (int i = 0; i < height; i++) {
            double[] newDist = new double[width];
            for (int j = 0; j < width; j++) {
                if (i == 0) {
                    distTo[j] = energyGrid[j][i];
                } else {
                    double left, middle, right;
                    if (j == 0) left = 1000*height;
                    else left = distTo[j - 1];
                    middle = distTo[j];
                    if (j == width - 1) right = 1000*height;
                    else right = distTo[j + 1];
                    edgeTo[j][i] = j + min(left, middle, right);
                    newDist[j] = energyGrid[j][i] + distTo[edgeTo[j][i]];
                }
            }
            distTo = newDist;
        }
        double minDist = 1000*height;
        int[] minRoute = new int[height];
        for (int j = 0; j < width; j++) {
            if (distTo[j] < minDist) {
                minDist = distTo[j];
                minRoute[height - 1] = j;
            }
        }

        for (int i = height - 2; i >= 0; i--) {
            minRoute[i] = edgeTo[minRoute[i + 1]][i + 1];
        }

        return minRoute;
    }

    private int min(double left, double middle, double right) {
        double min = middle;
        int result = 0;
        if (left < min) {
            min = left; 
            result = -1;
        }
        if (right < min) {
            result = 1;
        }
        return result;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (height <= 1) throw new IllegalArgumentException();
        if (seam.length != width) throw new IllegalArgumentException();
        height--;
        for (int i = 0; i < width; i++) {
            if (seam[i] < 0 || seam[i] > height) 
                throw new IllegalArgumentException("seam[" + i + "]: " + seam[i]);
            if (i > 0 && Math.abs(seam[i] - seam[i - 1]) > 1)
                throw new IllegalArgumentException();

            for (int j = seam[i]; j < height; j++) {
                pic.set(i, j, pic.get(i, j + 1));
                energyGrid[i][j] = energyGrid[i][j + 1];
            }
        }
        for (int i = 0; i < width; i++) {
            if (seam[i] < height) {
                energyGrid[i][seam[i]] = energyCal(i, seam[i]);
                if (i > 0) 
                    energyGrid[i - 1][seam[i]] = energyCal(i - 1, seam[i]);
                if (i < width - 1) 
                    energyGrid[i + 1][seam[i]] = energyCal(i + 1, seam[i]);
            }
            if (seam[i] > 0) 
                energyGrid[i][seam[i] - 1] = energyCal(i, seam[i] - 1);
        }
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (width <= 1) throw new IllegalArgumentException();
        if (seam.length != height) throw new IllegalArgumentException();
        width--;
        for (int i = 0; i < height; i++) {
            if (seam[i] < 0 || seam[i] > width) 
                throw new IllegalArgumentException("seam[" + i + "]: " + seam[i]);
            if (i > 0 && Math.abs(seam[i] - seam[i - 1]) > 1)
                throw new IllegalArgumentException();

            for (int j = seam[i]; j < width; j++) {
                pic.set(j, i, pic.get(j + 1, i));
                energyGrid[j][i] = energyGrid[j + 1][i];
            }
        }
        for (int i = 0; i < height; i++) {
            if (seam[i] < width) {
                energyGrid[seam[i]][i] = energyCal(seam[i], i);
                if (i > 0) 
                    energyGrid[seam[i]][i - 1] = energyCal(seam[i], i - 1);
                if (i < height - 1) 
                    energyGrid[seam[i]][i + 1] = energyCal(seam[i], i + 1);
            }
            if (seam[i] > 0) 
                energyGrid[seam[i] - 1][i] = energyCal(seam[i] - 1, i);
        }
    }
}