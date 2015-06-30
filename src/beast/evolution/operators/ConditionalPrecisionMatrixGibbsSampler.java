package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.math.distributions.Wishart;
import beast.math.matrixAlgebra.Matrix;

import java.lang.reflect.Array;

/**
 * @author Chieh-Hsi Wu
 */
public class ConditionalPrecisionMatrixGibbsSampler extends Operator {
    public Input<Integer> dfInput = new Input<Integer>("df", "The degrees of freedom of the Wishart distribution.", Input.Validate.REQUIRED);
    public Input<RealParameter> scaleMatrixInput = new Input<RealParameter>("scaleMatrix", "The scale matrix (as a vector) of the Wishart distribution.", Input.Validate.REQUIRED);
    public Input<ParameterList> xInput = new Input<ParameterList>("x","A list of values assumed to be from a multivariate normal distribution with the precision matrix to be sampled by this operator.");
    public Input<RealParameter> precisionInput = new Input<RealParameter>("precision","The precision matrix of a multivariate normal distribution.");
    public Input<RealParameter> meanInput = new Input<RealParameter>("mean", "The mean of a multivariate normal distribution.");
    public Input<ParameterList> conditionInput = new Input<ParameterList>(
            "condition",
            "The values that the distribution of x is conditioned on.",
            Input.Validate.REQUIRED
    );

    public Input<Double> correctConditionInput = new Input<Double>(
            "correctCondition",
            "The correct condition of the distribution.",
            Input.Validate.REQUIRED
    );
    private int df;
    private double[][] scaleMatrixInv;
    private ParameterList x;
    private int dim;
    private RealParameter mean;
    private ParameterList condition;
    private double correctCondition;

    public void initAndValidate(){
        df = dfInput.get();
        RealParameter scaleMatrixParam = scaleMatrixInput.get();
        dim = (int)(Math.sqrt(scaleMatrixParam.getDimension()));
        scaleMatrixInv = new double[dim][dim];
        int k = 0;
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                scaleMatrixInv[i][j] = scaleMatrixParam.getValue(k++);

            }

        }

        scaleMatrixInv = (new Matrix(scaleMatrixInv)).inverse().toComponents();
        mean = meanInput.get();
        x = xInput.get();
        condition = conditionInput.get();
        correctCondition = correctConditionInput.get();
    }

    public double proposal(){

        int xCount = x.getDimension();
        double[][] centredObs = new double[xCount][dim];
        for(int j = 0; j < dim; j++){
            double colMean = mean.getValue(j);
            for(int i = 0; i < xCount; i++){
                if(condition.getValue(i) == correctCondition){
                    centredObs[i][j] = x.getValue(i,j) - colMean;

                }

            }
        }

        int obsCount = 0;
        for(int  i = 0; i < xCount; i++){
            if(condition.getValue(i) == correctCondition){
                obsCount++;
            }
        }



        double[][] newScale = new double[dim][dim];
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                for(int k = 0; k < xCount; k++){
                    //System.out.println(i+" "+j+" "+k);
                    //if(condition.getValue(i) == correctCondition){
                    newScale[i][j] += centredObs[k][i]*centredObs[k][j];
                }

                newScale[i][j] +=scaleMatrixInv[i][j];
                //System.out.print(newScale[i][j]+" ");
            }
            //System.out.println();
        }
        Matrix newScaleMatrix = (new Matrix(newScale)).inverse();
        newScale = newScaleMatrix.toComponents();
        Double[] newScale1d = new Double[dim*dim];
        int u = 0;
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                newScale1d[u++] = newScale[i][j];
            }
        }


        double[][] draw = Wishart.nextWishart(df + obsCount, newScale);
        RealParameter precision = precisionInput.get(this);


        double logq = 0;
        try{
            Wishart w = new Wishart();
            w.initByName(
                    "scaleMatrix", new RealParameter(newScale1d),
                    "df", df + obsCount
            );

            logq = w.calcLogP(precision);
            //System.out.println("flag1: "+w.calcLogP(precision));
            int k = 0;
            for(int i = 0; i < dim; i ++){
                for(int j = 0; j < dim; j ++){
                //System.out.print(newScale[i][j]+" ");
                    precision.setValue(k++, draw[i][j]);
                }
            //System.out.println();

            }

            logq -= w.calcLogP(precision);
            //System.out.println("flag2: "+w.calcLogP(precision));
        }catch(Exception e){
            throw new RuntimeException(e);

        }



        return Double.POSITIVE_INFINITY;
        //return logq;

    }

    public static void main(String[] args){
        try{
            int df = 5;
            RealParameter scaleMatrix = new RealParameter(new Double[]{
                    5.11210762331838619588, -0.26905829596412556004, -1.25560538116591935420,
                    -0.26905829596412561555, 2.64573991031390098883, -0.98654708520179368314,
                    -1.25560538116591935420, -0.98654708520179357212, 2.06278026905829614535
            });
            RealParameter x1 = new RealParameter(new Double[]{1.0,3.0,1.0});
            RealParameter x2 = new RealParameter(new Double[]{2.0,1.0,5.0});
            RealParameter x3 = new RealParameter(new Double[]{2.0,4.0,2.0});
            ParameterList x = new ParameterList();
            x.initByName(
                    "parameter", x1,
                    "parameter", x2,
                    "parameter", x3
            );

            RealParameter precision = new RealParameter(
                    new Double[]{
                            0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0
                    }
            );
            RealParameter mean = new RealParameter(new Double[]{1.0, 2.0, 3.0});

            PrecisionWishartSampler pws = new PrecisionWishartSampler();
            pws.initByName(
                    "df", df,
                    "scaleMatrix", scaleMatrix,
                    "x", x,
                    "precision", precision,
                    "mean", mean
            );
            pws.proposal();
                /*
                 0.9278745409071454 -0.527703296455691 -0.4227034091165142
                 -0.5277032964556911 0.6158040602960733 0.4249566255830202
                 -0.4227034091165142 0.4249566255830202 0.4030253186423623
                 */
                /*Calculation by hand:
                0.9278745409071454, -0.5277032964556909, -0.4227034091165141,
                -0.5277032964556909,  0.6158040602960730,  0.4249566255830201,
                -0.4227034091165142,  0.4249566255830201,  0.4030253186423622 */

        }catch(Exception e){
            throw new RuntimeException(e);

        }


    }
}
