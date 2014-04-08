package beast.evolution.likelihood;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.GeneralUnitAlignment;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 17/07/13
 * Time: 12:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralUnitSepTempTreeLikelihood extends SepTempTreeLikelihood{

    public Input<ExtendedSepTempWVTreeLikelihood> sepTempWVTreeLikelihoodInput = new Input<ExtendedSepTempWVTreeLikelihood>(
            "sepTempTreeLikelihood",
            "The weight variable temporary tree likelihood for seperate subst model and rate category structure.",
            Input.Validate.REQUIRED
    );

    private ExtendedSepTempWVTreeLikelihood sepTempWVTreeLikelihood;
    public GeneralUnitSepTempTreeLikelihood(){
        siteModelInput.setRule(Input.Validate.OPTIONAL);
        treeInput.setRule(Input.Validate.OPTIONAL);
        dpNtdRateSepSiteModelInput.setRule(Input.Validate.OPTIONAL);
    }
    public void initAndValidate(){
        sepTempWVTreeLikelihood = sepTempWVTreeLikelihoodInput.get();
        alignment = dataInput.get();
    }
    public double calculateLogP(int unitIndex){
        throw new RuntimeException("Method to do");
    }

    public double calculateLogP(
            RealParameter rateParameter,
            int unitIndex){

        try{
            GeneralUnitAlignment guAlignment = (GeneralUnitAlignment)alignment;
            int[] sites = guAlignment.getSitesByUnit(unitIndex);
            double[] siteLogPs = sepTempWVTreeLikelihood.calculateLogP(rateParameter, sites);
            logP = 0.0;
            for(double siteLogP: siteLogPs){
                logP+=siteLogP;
            }
            return logP;
        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }


    public double calculateLogP(
            double alpha,
            double invPr,
            double rate,
            double siteModelChoice,
            int unitIndex){

        try{
            GeneralUnitAlignment guAlignment = (GeneralUnitAlignment)alignment;
            int[] sites = guAlignment.getSitesByUnit(unitIndex);
            double[] siteLogPs = sepTempWVTreeLikelihood.calculateLogP(
                    alpha,
                    invPr,
                    rate,
                    siteModelChoice,
                    sites);
            logP = 0.0;
            for(double siteLogP: siteLogPs){
                logP+=siteLogP;
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
            int unitIndex){
        try{

            GeneralUnitAlignment guAlignment = (GeneralUnitAlignment)alignment;
            int[] sites = guAlignment.getSitesByUnit(unitIndex);
            double[] siteLogPs = sepTempWVTreeLikelihood.calculateLogP(
                    modelParameters,
                    modelCode,
                    freqs,
                    sites);
            logP = 0.0;
            for(double siteLogP: siteLogPs){
                logP+=siteLogP;
            }
            return logP;


        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }




}
