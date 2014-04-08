package beast.evolution.likelihood;

import beast.core.parameter.RealParameter;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.sitemodel.QuietGammaSiteBMA;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

/**
 * @author Chieh-Hsi Wu
 */
public class ExtendedTempWVTreeLikelihood extends TempWVTreeLikelihood {
    public double[] calculateLogP (
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter alpha,
            RealParameter invPr,
            RealParameter rate,
            RealParameter siteModelIndicator,
            int[] sites) throws Exception{
        double[] siteLogP = new double[sites.length];
        calculateLogP(modelParameters,modelCode,freqs,alpha, invPr,rate, siteModelIndicator);
        for(int i = 0; i < sites.length;i++){
            siteLogP[i] = m_fPatternLogLikelihoods[m_data.get().getPatternIndex(sites[i])];
        }
        return siteLogP;
    }

    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter alpha,
            RealParameter invPr,
            RealParameter rate,
            RealParameter siteModelIndicator) throws Exception{


        SwitchingNtdBMA substModel = (SwitchingNtdBMA)m_substitutionModel;
        (substModel.getLogKappa()).setValueQuietly(0, modelParameters.getValue(0));
        (substModel.getLogTN()).setValueQuietly(0,modelParameters.getValue(1));
        (substModel.getLogAC()).setValueQuietly(0,modelParameters.getValue(2));
        (substModel.getLogAT()).setValueQuietly(0,modelParameters.getValue(3));
        (substModel.getLogGC()).setValueQuietly(0,modelParameters.getValue(4));
        (substModel.getModelChoose()).setValueQuietly(0,modelCode.getValue());
        (substModel.getFreqs()).setValueQuietly(0,freqs.getValue(0));
        (substModel.getFreqs()).setValueQuietly(1,freqs.getValue(1));
        (substModel.getFreqs()).setValueQuietly(2,freqs.getValue(2));
        (substModel.getFreqs()).setValueQuietly(3,freqs.getValue(3));
        substModel.setUpdateMatrix(true);
        //((QuietGammaSiteBMA)m_siteModel).getRateParameter().setValueQuietly(0,rate.getValue());
        ((QuietGammaSiteBMA)m_siteModel).setMuValueQuietly(rate.getValue());
        ((QuietGammaSiteBMA)m_siteModel).setShapeValueQuietly(alpha.getValue());
        ((QuietGammaSiteBMA)m_siteModel).setInvPrValueQuietly(invPr.getValue());
        ((QuietGammaSiteBMA)m_siteModel).setModelChoiceQuietly(siteModelIndicator.getValue());
        ((QuietGammaSiteBMA)m_siteModel).setRatesKnown(false);
        return calculateLogP();



    }



}
