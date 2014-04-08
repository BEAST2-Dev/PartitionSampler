package beast.evolution.likelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.substitutionmodel.NtdBMA;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 15/07/13
 * Time: 6:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralUnitTempTreeLikelihood extends  TempTreeLikelihood{
    public Input<TempWVTreeLikelihood> tempWVTreeLikelihoodInput = new Input<TempWVTreeLikelihood>(
            "tempWVTreeLikelihood",
            "The weight variable temporary tree likelihood.",
            Input.Validate.REQUIRED
    );

    public GeneralUnitTempTreeLikelihood(){
        siteModelInput.setRule(Input.Validate.OPTIONAL);
        treeInput.setRule(Input.Validate.OPTIONAL);
    }

    private TempWVTreeLikelihood tempWVTreeLikelihood;
    public void initAndValidate(){
        alignment = dataInput.get();
        tempWVTreeLikelihood = tempWVTreeLikelihoodInput.get();
    }
    public double calculateLogP(int unitIndex){
        throw new RuntimeException("Method to do");
    }

    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int unitIndex){
        try{

            if(substModel instanceof NtdBMA){
                ((NtdBMA)substModel).getLogKappa().setValueQuietly(0,modelParameters.getValue(0));
                ((NtdBMA)substModel).getLogTN().setValueQuietly(0,modelParameters.getValue(1));
                ((NtdBMA)substModel).getLogAC().setValueQuietly(0,modelParameters.getValue(2));
                ((NtdBMA)substModel).getLogAT().setValueQuietly(0,modelParameters.getValue(3));
                ((NtdBMA)substModel).getLogGC().setValueQuietly(0,modelParameters.getValue(4));
                ((NtdBMA)substModel).getModelChoose().setValueQuietly(0,modelCode.getValue());
                ((NtdBMA)substModel).getFreqs().setValueQuietly(0,freqs.getValue(0));
                ((NtdBMA)substModel).getFreqs().setValueQuietly(1,freqs.getValue(1));
                ((NtdBMA)substModel).getFreqs().setValueQuietly(2,freqs.getValue(2));
                ((NtdBMA)substModel).getFreqs().setValueQuietly(3,freqs.getValue(3));

            }else{
                throw new RuntimeException("Need NtdBMA");
            }

            ((NtdBMA)substModel).setUpdateMatrix(true);

            int[] siteindices = ((GeneralUnitAlignment)alignment).getSitesByUnit(unitIndex);
            tempWVTreeLikelihood.setupPatternWeightsFromSites(siteindices);
            double[] siteLogPs = tempWVTreeLikelihood.calculateLogP(modelParameters,modelCode,freqs,siteindices);
            logP = 0;
            for(double siteLogP: siteLogPs){
                logP += siteLogP;
            }
            return logP;
        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }



    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter rate,
            int unitIndex){
        try{

            //System.out.println("Hi");

            int[] siteindices = ((GeneralUnitAlignment)alignment).getSitesByUnit(unitIndex);
            tempWVTreeLikelihood.setupPatternWeightsFromSites(siteindices);
            double[] siteLogPs = tempWVTreeLikelihood.calculateLogP(modelParameters,modelCode,freqs, rate,siteindices);
            logP = 0;
            for(double siteLogP: siteLogPs){
                //System.out.println(siteLogP);
                logP += siteLogP;
            }

            //System.out.println("logP: "+logP);
            return logP;


        }catch(Exception e){
            throw new RuntimeException(e);

        }
    }


    public double calculateLogP(
            RealParameter rateParameter,
            int unitIndex){
        try{
            siteModelInput.get().getRateParameter().setValueQuietly(0,rateParameter.getValue());


            int[] siteindices = ((GeneralUnitAlignment)alignment).getSitesByUnit(unitIndex);
            tempWVTreeLikelihood.setupPatternWeightsFromSites(siteindices);
            double[] siteLogPs = tempWVTreeLikelihood.calculateLogP(rateParameter,siteindices);
            logP = 0;
            for(double siteLogP: siteLogPs){
                logP += siteLogP;
            }
            return logP;
        }catch(Exception e){
            throw new RuntimeException(e);

        }
    }


    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter alpha,
            RealParameter invPr,
            RealParameter rate,
            RealParameter siteModelIndicator,
            int unitIndex){
        try{


            int[] siteindices = ((GeneralUnitAlignment)alignment).getSitesByUnit(unitIndex);
            tempWVTreeLikelihood.setupPatternWeightsFromSites(siteindices);
            double[] siteLogPs = ((ExtendedTempWVTreeLikelihood)tempWVTreeLikelihood).calculateLogP(
                    modelParameters,
                    modelCode,
                    freqs,
                    alpha,
                    invPr,
                    rate,
                    siteModelIndicator,
                    siteindices);
            logP = 0;
            for(double siteLogP: siteLogPs){
                //System.out.println(siteLogP);
                logP += siteLogP;
            }

            //System.out.println("logP: "+logP);
            return logP;


        }catch(Exception e){
            throw new RuntimeException(e);

        }
    }



}
