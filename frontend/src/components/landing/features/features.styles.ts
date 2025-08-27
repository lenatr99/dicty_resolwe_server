import styled from 'styled-components';
import { breakpoints } from 'components/app/globalStyle';

export const FeaturesWrapper = styled.div`
    width: 100%;
    display: flex;
    flex-flow: column nowrap;
    padding-bottom: 35px;

    @media (min-width: ${breakpoints.small}px) {
        flex-flow: row wrap;
        justify-content: space-between;
    }
`;

export const FeatureWrapper = styled.div`
    width: 100%;
    text-align: center;
    display: flex;
    flex-direction: column;
    align-content: center;

    @media (min-width: ${breakpoints.small}px) {
        width: 24%;
    }
`;

export const FeatureDescription = styled.p`
    padding-top: 20px;
`;

export const FeatureImage = styled.img`
    object-fit: none;
    max-width: 100%;
`;
