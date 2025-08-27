import styled, { css } from 'styled-components';
import { breakpoints } from 'components/app/globalStyle';

type ContainerProps = {
    $centerText?: boolean;
    $paddingTop?: number;
    $paddingBottom?: number;
};

export const DarkSectionContainer = styled.div`
    background-color: #f5f5f5;
    width: 100%;
`;

export const ContentContainer = styled.div`
    margin-right: auto;
    margin-left: auto;
    height: 100%;
    width: 100%;

    @media (min-width: ${breakpoints.small}px) {
        width: 750px;
    }
    @media (min-width: ${breakpoints.mid}px) {
        width: 970px;
    }
    @media (min-width: ${breakpoints.large}px) {
        width: 1170px;
    }
`;

export const SectionContentContainer = styled(ContentContainer)<ContainerProps>`
    padding-left: 15px;
    padding-right: 15px;
    padding-top: ${(props): string => `${props.$paddingTop ?? 35}px`};
    padding-bottom: ${(props): string => `${props.$paddingBottom ?? 35}px`};

    ${(props) =>
        props.$centerText
            ? css`
                  text-align: center;
              `
            : null};
`;

export const SectionHorizontalLine = styled.hr`
    color: #6f6f6f;
    border: 0;
    border-top: 1px solid #eee;
    margin-right: auto;
    margin-left: auto;

    @media (min-width: ${breakpoints.small}px) {
        width: 750px;
    }
    @media (min-width: ${breakpoints.mid}px) {
        width: 970px;
    }
    @media (min-width: ${breakpoints.large}px) {
        width: 1170px;
    }
`;

export const AlignCenter = styled.div`
    display: flex;
    flex-flow: row nowrap;
    justify-content: center;
`;
